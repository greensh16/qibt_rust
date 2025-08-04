use crate::config::StreamingConfig;
use crate::io::cloud_auth::{CloudAuthProvider, CloudAuthError};
use std::collections::HashMap;
use std::sync::Arc;

#[cfg(feature = "zarr")]
use {
    bytes::Bytes,
    tokio::sync::RwLock,
    crate::io::cloud_auth::{RetryConfig, retry_with_backoff},
};

#[cfg(feature = "zarr")]
use object_store::{ObjectStore, path::Path as ObjectPath, GetOptions, GetRange};

/// Streaming reader with range request support and chunk-wise data access
pub struct StreamingReader {
    #[cfg(feature = "zarr")]
    store: Option<Arc<dyn ObjectStore>>,
    config: StreamingConfig,
    auth_provider: CloudAuthProvider,
    #[cfg(feature = "zarr")]
    cache: Arc<RwLock<HashMap<String, CachedChunk>>>,
}

/// Cached chunk of data
#[cfg(feature = "zarr")]
#[derive(Clone)]
struct CachedChunk {
    data: Bytes,
    start_offset: u64,
    end_offset: u64,
    last_accessed: std::time::Instant,
}

/// Range specification for data requests
#[derive(Debug, Clone)]
pub struct DataRange {
    pub start: u64,
    pub end: Option<u64>, // None means read to end
}

impl DataRange {
    pub fn new(start: u64, end: Option<u64>) -> Self {
        Self { start, end }
    }
    
    pub fn from_offset_and_length(offset: u64, length: u64) -> Self {
        Self {
            start: offset,
            end: Some(offset + length),
        }
    }
    
    pub fn length(&self) -> Option<u64> {
        self.end.map(|end| end - self.start)
    }
}

/// Error types for streaming operations
#[derive(Debug, thiserror::Error)]
pub enum StreamingError {
    #[error("Cloud authentication error: {0}")]
    CloudAuth(#[from] CloudAuthError),
    
    #[cfg(feature = "zarr")]
    #[error("Object store error: {0}")]
    ObjectStore(#[from] object_store::Error),
    
    #[error("Invalid range: {0}")]
    InvalidRange(String),
    
    #[error("Cache error: {0}")]
    CacheError(String),
    
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

impl StreamingReader {
    /// Create a new streaming reader
    pub fn new(
        config: StreamingConfig,
        auth_provider: CloudAuthProvider,
    ) -> Self {
        Self {
            #[cfg(feature = "zarr")]
            store: None,
            config,
            auth_provider,
            #[cfg(feature = "zarr")]
            cache: Arc::new(RwLock::new(HashMap::new())),
        }
    }

    /// Initialize the reader with a remote URL
    #[cfg(feature = "zarr")]
    pub async fn initialize(&mut self, url: &str) -> Result<(), StreamingError> {
        let store = self.auth_provider.create_store_from_url(url)?;
        self.store = Some(Arc::new(store));
        Ok(())
    }

    #[cfg(not(feature = "zarr"))]
    pub async fn initialize(&mut self, _url: &str) -> Result<(), StreamingError> {
        Err(StreamingError::CacheError(
            "Zarr feature not enabled".to_string()
        ))
    }

    /// Read data from a remote object with range support
    #[cfg(feature = "zarr")]
    pub async fn read_range(
        &self,
        object_path: &str,
        range: DataRange,
    ) -> Result<Bytes, StreamingError> {
        let store = self.store.as_ref()
            .ok_or_else(|| StreamingError::CacheError("Store not initialized".to_string()))?;

        // Check cache first
        if let Some(cached_data) = self.check_cache(object_path, &range).await {
            return Ok(cached_data);
        }

        let path = ObjectPath::from(object_path);
        let retry_config = RetryConfig::from(&self.config);

        let get_options = if self.config.enable_range_requests {
            let range_spec = match range.end {
                Some(end) => GetRange::Bounded(range.start as usize, end as usize),
                None => GetRange::Offset(range.start as usize),
            };
            GetOptions {
                range: Some(range_spec),
                ..Default::default()
            }
        } else {
            GetOptions::default()
        };

        // Retry logic for network operations
        let data = retry_with_backoff(
            || async {
                store.get_opts(&path, get_options.clone()).await
            },
            &retry_config,
        ).await?;

        let bytes = data.bytes().await?;

        // Cache the data
        self.cache_data(object_path, &range, bytes.clone()).await;

        Ok(bytes)
    }

    #[cfg(not(feature = "zarr"))]
    pub async fn read_range(
        &self,
        _object_path: &str,
        _range: DataRange,
    ) -> Result<Bytes, StreamingError> {
        Err(StreamingError::CacheError(
            "Zarr feature not enabled".to_string()
        ))
    }

    /// Read data in chunks for streaming processing
    #[cfg(feature = "zarr")]
    pub async fn read_chunked<F>(
        &self,
        object_path: &str,
        chunk_processor: F,
    ) -> Result<(), StreamingError>
    where
        F: Fn(Bytes, u64) -> Result<(), StreamingError>,
    {
        let store = self.store.as_ref()
            .ok_or_else(|| StreamingError::CacheError("Store not initialized".to_string()))?;

        let path = ObjectPath::from(object_path);
        
        // Get object metadata to determine size
        let metadata = store.head(&path).await?;
        let total_size = metadata.size;
        
        let chunk_size = self.config.chunk_size as u64;
        let mut offset = 0;

        while offset < total_size {
            let end = std::cmp::min(offset + chunk_size, total_size);
            let range = DataRange::new(offset, Some(end));
            
            let chunk_data = self.read_range(object_path, range).await?;
            chunk_processor(chunk_data, offset)?;
            
            offset = end;
        }

        Ok(())
    }

    #[cfg(not(feature = "zarr"))]
    pub async fn read_chunked<F>(
        &self,
        _object_path: &str,
        _chunk_processor: F,
    ) -> Result<(), StreamingError>
    where
        F: Fn(Bytes, u64) -> Result<(), StreamingError>,
    {
        Err(StreamingError::CacheError(
            "Zarr feature not enabled".to_string()
        ))
    }

    /// Read multiple ranges concurrently
    #[cfg(feature = "zarr")]
    pub async fn read_multiple_ranges(
        &self,
        object_path: &str,
        ranges: Vec<DataRange>,
    ) -> Result<Vec<Bytes>, StreamingError> {
        use futures::future::try_join_all;
        
        // Limit concurrency to avoid overwhelming the server
        let semaphore = Arc::new(tokio::sync::Semaphore::new(
            self.config.max_concurrent_requests
        ));
        
        let futures = ranges.into_iter().map(|range| {
            let semaphore = semaphore.clone();
            let object_path = object_path.to_string();
            async move {
                let _permit = semaphore.acquire().await.unwrap();
                self.read_range(&object_path, range).await
            }
        });

        try_join_all(futures).await
    }

    #[cfg(not(feature = "zarr"))]
    pub async fn read_multiple_ranges(
        &self,
        _object_path: &str,
        _ranges: Vec<DataRange>,
    ) -> Result<Vec<Bytes>, StreamingError> {
        Err(StreamingError::CacheError(
            "Zarr feature not enabled".to_string()
        ))
    }

    /// Get object metadata (size, last modified, etc.)
    #[cfg(feature = "zarr")]
    pub async fn get_object_metadata(
        &self,
        object_path: &str,
    ) -> Result<object_store::ObjectMeta, StreamingError> {
        let store = self.store.as_ref()
            .ok_or_else(|| StreamingError::CacheError("Store not initialized".to_string()))?;

        let path = ObjectPath::from(object_path);
        let metadata = store.head(&path).await?;
        Ok(metadata)
    }

    #[cfg(not(feature = "zarr"))]
    pub async fn get_object_metadata(
        &self,
        _object_path: &str,
    ) -> Result<(), StreamingError> {
        Err(StreamingError::CacheError(
            "Zarr feature not enabled".to_string()
        ))
    }

    /// Check if data is available in cache
    #[cfg(feature = "zarr")]
    async fn check_cache(&self, object_path: &str, range: &DataRange) -> Option<Bytes> {
        let cache = self.cache.read().await;
        let cache_key = format!("{}:{}:{}", object_path, range.start, range.end.unwrap_or(u64::MAX));
        
        if let Some(cached_chunk) = cache.get(&cache_key) {
            // Check if cached chunk covers the requested range
            if cached_chunk.start_offset <= range.start {
                let end_offset = range.end.unwrap_or(cached_chunk.end_offset);
                if cached_chunk.end_offset >= end_offset {
                    let start_pos = (range.start - cached_chunk.start_offset) as usize;
                    let end_pos = (end_offset - cached_chunk.start_offset) as usize;
                    return Some(cached_chunk.data.slice(start_pos..end_pos));
                }
            }
        }
        
        None
    }

    /// Cache data for future requests
    #[cfg(feature = "zarr")]
    async fn cache_data(&self, object_path: &str, range: &DataRange, data: Bytes) {
        if data.len() > self.config.buffer_size {
            // Don't cache chunks larger than buffer size
            return;
        }

        let cache_key = format!("{}:{}:{}", object_path, range.start, range.end.unwrap_or(u64::MAX));
        let end_offset = range.end.unwrap_or(range.start + data.len() as u64);
        
        let cached_chunk = CachedChunk {
            data,
            start_offset: range.start,
            end_offset,
            last_accessed: std::time::Instant::now(),
        };

        let mut cache = self.cache.write().await;
        
        // Simple cache eviction: remove oldest entries if cache is full
        if cache.len() >= 1000 { // Max 1000 cached chunks
            let oldest_key = cache.iter()
                .min_by_key(|(_, chunk)| chunk.last_accessed)
                .map(|(key, _)| key.clone());
            
            if let Some(key) = oldest_key {
                cache.remove(&key);
            }
        }
        
        cache.insert(cache_key, cached_chunk);
    }

    /// Clear the cache
    #[cfg(feature = "zarr")]
    pub async fn clear_cache(&self) {
        let mut cache = self.cache.write().await;
        cache.clear();
    }
    
    #[cfg(not(feature = "zarr"))]
    pub async fn clear_cache(&self) {
        // No-op when zarr feature is disabled
    }

    /// Get cache statistics
    #[cfg(feature = "zarr")]
    pub async fn get_cache_stats(&self) -> (usize, usize) {
        let cache = self.cache.read().await;
        let num_entries = cache.len();
        let total_size: usize = cache.values()
            .map(|chunk| chunk.data.len())
            .sum();
        (num_entries, total_size)
    }
    
    #[cfg(not(feature = "zarr"))]
    pub async fn get_cache_stats(&self) -> (usize, usize) {
        (0, 0)
    }
}

/// Helper trait for converting byte ranges to data types
pub trait FromBytes {
    fn from_bytes(data: &[u8]) -> Result<Self, StreamingError>
    where
        Self: Sized;
}

impl FromBytes for Vec<f32> {
    fn from_bytes(data: &[u8]) -> Result<Self, StreamingError> {
        if data.len() % 4 != 0 {
            return Err(StreamingError::InvalidRange(
                "Data length not aligned to f32 boundary".to_string()
            ));
        }
        
        let mut result = Vec::with_capacity(data.len() / 4);
        for chunk in data.chunks_exact(4) {
            let bytes: [u8; 4] = chunk.try_into().unwrap();
            result.push(f32::from_le_bytes(bytes));
        }
        
        Ok(result)
    }
}

impl FromBytes for Vec<f64> {
    fn from_bytes(data: &[u8]) -> Result<Self, StreamingError> {
        if data.len() % 8 != 0 {
            return Err(StreamingError::InvalidRange(
                "Data length not aligned to f64 boundary".to_string()
            ));
        }
        
        let mut result = Vec::with_capacity(data.len() / 8);
        for chunk in data.chunks_exact(8) {
            let bytes: [u8; 8] = chunk.try_into().unwrap();
            result.push(f64::from_le_bytes(bytes));
        }
        
        Ok(result)
    }
}

#[cfg(all(test, feature = "zarr"))]
mod tests {
    use super::*;
    use crate::config::{CloudAuth, StreamingConfig};

    #[test]
    fn test_data_range() {
        let range1 = DataRange::new(0, Some(100));
        assert_eq!(range1.length(), Some(100));
        
        let range2 = DataRange::new(50, None);
        assert_eq!(range2.length(), None);
        
        let range3 = DataRange::from_offset_and_length(10, 50);
        assert_eq!(range3.start, 10);
        assert_eq!(range3.end, Some(60));
        assert_eq!(range3.length(), Some(50));
    }

    #[test]
    fn test_from_bytes_f32() {
        let data = vec![0u8, 0, 0x80, 0x3F]; // 1.0 in little-endian f32
        let result: Vec<f32> = FromBytes::from_bytes(&data).unwrap();
        assert_eq!(result.len(), 1);
        assert!((result[0] - 1.0).abs() < f32::EPSILON);
    }

    #[test]
    fn test_from_bytes_f64() {
        let data = vec![0u8, 0, 0, 0, 0, 0, 0xF0, 0x3F]; // 1.0 in little-endian f64
        let result: Vec<f64> = FromBytes::from_bytes(&data).unwrap();
        assert_eq!(result.len(), 1);
        assert!((result[0] - 1.0).abs() < f64::EPSILON);
    }

    #[tokio::test]
    async fn test_streaming_reader_creation() {
        let config = StreamingConfig::default();
        let auth_provider = CloudAuthProvider::new(
            CloudAuth::default(),
            config.clone(),
        );
        
        let reader = StreamingReader::new(config, auth_provider);
        let (entries, size) = reader.get_cache_stats().await;
        assert_eq!(entries, 0);
        assert_eq!(size, 0);
    }

    #[tokio::test]
    async fn test_cache_operations() {
        let config = StreamingConfig::default();
        let auth_provider = CloudAuthProvider::new(
            CloudAuth::default(),
            config.clone(),
        );
        
        let reader = StreamingReader::new(config, auth_provider);
        
        // Cache some data
        let test_data = Bytes::from("test data");
        let range = DataRange::new(0, Some(9));
        reader.cache_data("test/path", &range, test_data.clone()).await;
        
        let (entries, size) = reader.get_cache_stats().await;
        assert_eq!(entries, 1);
        assert_eq!(size, 9);
        
        // Check cache hit
        let cached = reader.check_cache("test/path", &range).await;
        assert!(cached.is_some());
        assert_eq!(cached.unwrap(), test_data);
        
        // Clear cache
        reader.clear_cache().await;
        let (entries, size) = reader.get_cache_stats().await;
        assert_eq!(entries, 0);
        assert_eq!(size, 0);
    }
}
