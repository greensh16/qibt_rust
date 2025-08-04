use crate::config::{CloudAuth, StreamingConfig};
use std::collections::HashMap;
use std::env;
use std::time::Duration;

#[cfg(feature = "zarr")]
use object_store::{
    ObjectStore, Result as ObjectStoreResult,
    aws::{AmazonS3Builder, AmazonS3},
    gcp::{GoogleCloudStorageBuilder, GoogleCloudStorage},
    azure::{MicrosoftAzureBuilder, MicrosoftAzure},
    http::{HttpBuilder, HttpStore},
};

/// Cloud storage error type
#[derive(Debug, thiserror::Error)]
pub enum CloudAuthError {
    #[error("Missing required environment variable: {0}")]
    MissingEnvVar(String),
    
    #[error("Invalid configuration: {0}")]
    InvalidConfig(String),
    
    #[error("Authentication failed: {0}")]
    AuthFailed(String),
    
    #[cfg(feature = "zarr")]
    #[error("Object store error: {0}")]
    ObjectStore(#[from] object_store::Error),
    
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

/// Environment variable loader with optional CLI/config override
pub struct CloudAuthProvider {
    auth_config: CloudAuth,
    streaming_config: StreamingConfig,
}

impl CloudAuthProvider {
    /// Create a new cloud authentication provider
    pub fn new(auth_config: CloudAuth, streaming_config: StreamingConfig) -> Self {
        Self {
            auth_config,
            streaming_config,
        }
    }

    /// Load authentication from environment variables or configuration
    pub fn load_auth(&self) -> Result<CloudAuth, CloudAuthError> {
        let mut auth = self.auth_config.clone();
        
        // AWS credentials - check env vars if not provided in config
        if auth.aws_access_key_id.is_none() {
            auth.aws_access_key_id = env::var("AWS_ACCESS_KEY_ID").ok()
                .or_else(|| env::var("QIBT_AWS_ACCESS_KEY_ID").ok());
        }
        
        if auth.aws_secret_access_key.is_none() {
            auth.aws_secret_access_key = env::var("AWS_SECRET_ACCESS_KEY").ok()
                .or_else(|| env::var("QIBT_AWS_SECRET_ACCESS_KEY").ok());
        }
        
        if auth.aws_session_token.is_none() {
            auth.aws_session_token = env::var("AWS_SESSION_TOKEN").ok()
                .or_else(|| env::var("QIBT_AWS_SESSION_TOKEN").ok());
        }
        
        if auth.aws_region.is_none() {
            auth.aws_region = env::var("AWS_REGION").ok()
                .or_else(|| env::var("AWS_DEFAULT_REGION").ok())
                .or_else(|| env::var("QIBT_AWS_REGION").ok());
        }
        
        // GCP credentials
        if auth.gcp_service_account_key.is_none() {
            auth.gcp_service_account_key = env::var("GOOGLE_APPLICATION_CREDENTIALS").ok()
                .or_else(|| env::var("QIBT_GCP_SERVICE_ACCOUNT_KEY").ok());
        }
        
        if auth.gcp_project_id.is_none() {
            auth.gcp_project_id = env::var("GOOGLE_CLOUD_PROJECT").ok()
                .or_else(|| env::var("QIBT_GCP_PROJECT_ID").ok());
        }
        
        // Azure credentials
        if auth.azure_storage_account.is_none() {
            auth.azure_storage_account = env::var("AZURE_STORAGE_ACCOUNT").ok()
                .or_else(|| env::var("QIBT_AZURE_STORAGE_ACCOUNT").ok());
        }
        
        if auth.azure_storage_key.is_none() {
            auth.azure_storage_key = env::var("AZURE_STORAGE_KEY").ok()
                .or_else(|| env::var("AZURE_STORAGE_ACCESS_KEY").ok())
                .or_else(|| env::var("QIBT_AZURE_STORAGE_KEY").ok());
        }
        
        // Custom endpoint
        if auth.endpoint_url.is_none() {
            auth.endpoint_url = env::var("QIBT_ENDPOINT_URL").ok()
                .or_else(|| env::var("AWS_ENDPOINT_URL").ok());
        }
        
        Ok(auth)
    }

    #[cfg(feature = "zarr")]
    /// Create object store for AWS S3
    pub fn create_s3_store(&self, bucket: &str, region: Option<&str>) -> Result<AmazonS3, CloudAuthError> {
        let auth = self.load_auth()?;
        
        let mut builder = AmazonS3Builder::new()
            .with_bucket_name(bucket);
        
        if let Some(access_key) = &auth.aws_access_key_id {
            builder = builder.with_access_key_id(access_key);
        }
        
        if let Some(secret_key) = &auth.aws_secret_access_key {
            builder = builder.with_secret_access_key(secret_key);
        }
        
        if let Some(session_token) = &auth.aws_session_token {
            builder = builder.with_token(session_token);
        }
        
        if let Some(region) = region.or(auth.aws_region.as_deref()) {
            builder = builder.with_region(region);
        }
        
        if let Some(endpoint) = &auth.endpoint_url {
            builder = builder.with_endpoint(endpoint);
        }
        
        // Configure client with retry and timeout settings
        let client_options = object_store::ClientOptions::new()
            .with_timeout(Duration::from_secs(self.streaming_config.request_timeout_secs))
            .with_connect_timeout(Duration::from_secs(10));
        
        builder = builder.with_client_options(client_options);
        
        builder.build().map_err(CloudAuthError::ObjectStore)
    }

    #[cfg(feature = "zarr")]
    /// Create object store for Google Cloud Storage
    pub fn create_gcs_store(&self, bucket: &str) -> Result<GoogleCloudStorage, CloudAuthError> {
        let auth = self.load_auth()?;
        
        let mut builder = GoogleCloudStorageBuilder::new()
            .with_bucket_name(bucket);
        
        if let Some(service_account_path) = &auth.gcp_service_account_key {
            builder = builder.with_service_account_path(service_account_path);
        }
        
        if let Some(project_id) = &auth.gcp_project_id {
            builder = builder.with_project_id(project_id);
        }
        
        // Configure client with retry and timeout settings
        let client_options = object_store::ClientOptions::new()
            .with_timeout(Duration::from_secs(self.streaming_config.request_timeout_secs))
            .with_connect_timeout(Duration::from_secs(10));
        
        builder = builder.with_client_options(client_options);
        
        builder.build().map_err(CloudAuthError::ObjectStore)
    }

    #[cfg(feature = "zarr")]
    /// Create object store for Azure Blob Storage
    pub fn create_azure_store(&self, account: &str, container: &str) -> Result<MicrosoftAzure, CloudAuthError> {
        let auth = self.load_auth()?;
        
        let mut builder = MicrosoftAzureBuilder::new()
            .with_account(account)
            .with_container_name(container);
        
        if let Some(access_key) = &auth.azure_storage_key {
            builder = builder.with_access_key(access_key);
        }
        
        // Configure client with retry and timeout settings
        let client_options = object_store::ClientOptions::new()
            .with_timeout(Duration::from_secs(self.streaming_config.request_timeout_secs))
            .with_connect_timeout(Duration::from_secs(10));
        
        builder = builder.with_client_options(client_options);
        
        builder.build().map_err(CloudAuthError::ObjectStore)
    }

    #[cfg(feature = "zarr")]
    /// Create HTTP object store with authentication headers
    pub fn create_http_store(&self, url: &str) -> Result<HttpStore, CloudAuthError> {
        let auth = self.load_auth()?;
        
        let mut builder = HttpBuilder::new()
            .with_url(url);
        
        // Add custom headers if provided
        for (key, value) in &auth.custom_headers {
            builder = builder.with_header(key, value);
        }
        
        // Configure client with retry and timeout settings
        let client_options = object_store::ClientOptions::new()
            .with_timeout(Duration::from_secs(self.streaming_config.request_timeout_secs))
            .with_connect_timeout(Duration::from_secs(10));
        
        builder = builder.with_client_options(client_options);
        
        builder.build().map_err(CloudAuthError::ObjectStore)
    }

    #[cfg(feature = "zarr")]
    /// Create object store based on URL scheme
    pub fn create_store_from_url(&self, url: &str) -> Result<Box<dyn ObjectStore>, CloudAuthError> {
        if url.starts_with("s3://") {
            let parts: Vec<&str> = url.strip_prefix("s3://").unwrap().splitn(2, '/').collect();
            if parts.len() < 1 {
                return Err(CloudAuthError::InvalidConfig("Invalid S3 URL format".to_string()));
            }
            let bucket = parts[0];
            let store = self.create_s3_store(bucket, None)?;
            Ok(Box::new(store))
        } else if url.starts_with("gs://") || url.starts_with("gcs://") {
            let parts: Vec<&str> = url.strip_prefix("gs://")
                .or_else(|| url.strip_prefix("gcs://"))
                .unwrap().splitn(2, '/').collect();
            if parts.len() < 1 {
                return Err(CloudAuthError::InvalidConfig("Invalid GCS URL format".to_string()));
            }
            let bucket = parts[0];
            let store = self.create_gcs_store(bucket)?;
            Ok(Box::new(store))
        } else if url.starts_with("abfs://") || url.starts_with("az://") {
            let parts: Vec<&str> = url.strip_prefix("abfs://")
                .or_else(|| url.strip_prefix("az://"))
                .unwrap().splitn(2, '/').collect();
            if parts.len() < 2 {
                return Err(CloudAuthError::InvalidConfig("Invalid Azure URL format".to_string()));
            }
            let container = parts[0];
            let account = parts[1].split('.').next().unwrap_or("");
            let store = self.create_azure_store(account, container)?;
            Ok(Box::new(store))
        } else if url.starts_with("http://") || url.starts_with("https://") {
            let store = self.create_http_store(url)?;
            Ok(Box::new(store))
        } else {
            Err(CloudAuthError::InvalidConfig(format!("Unsupported URL scheme: {}", url)))
        }
    }
}

/// Retry logic with exponential backoff
pub struct RetryConfig {
    pub max_retries: usize,
    pub base_delay_ms: u64,
    pub exponential_backoff: bool,
    pub max_delay_ms: u64,
}

impl From<&StreamingConfig> for RetryConfig {
    fn from(config: &StreamingConfig) -> Self {
        Self {
            max_retries: config.max_retries,
            base_delay_ms: config.retry_delay_ms,
            exponential_backoff: config.exponential_backoff,
            max_delay_ms: 30_000, // 30 seconds max delay
        }
    }
}

impl RetryConfig {
    /// Calculate delay for a given retry attempt
    pub fn calculate_delay(&self, attempt: usize) -> Duration {
        if attempt == 0 {
            return Duration::from_millis(0);
        }
        
        let delay_ms = if self.exponential_backoff {
            let exponential_delay = self.base_delay_ms * (2_u64.pow((attempt - 1) as u32));
            std::cmp::min(exponential_delay, self.max_delay_ms)
        } else {
            self.base_delay_ms
        };
        
        Duration::from_millis(delay_ms)
    }
}

/// Retry a fallible operation with exponential backoff
pub async fn retry_with_backoff<F, T, E>(
    mut operation: F,
    retry_config: &RetryConfig,
) -> Result<T, E>
where
    F: FnMut() -> Result<T, E>,
    E: std::fmt::Debug,
{
    let mut last_error = None;
    
    for attempt in 0..=retry_config.max_retries {
        match operation() {
            Ok(result) => return Ok(result),
            Err(error) => {
                last_error = Some(error);
                
                if attempt < retry_config.max_retries {
                    let delay = retry_config.calculate_delay(attempt + 1);
                    tokio::time::sleep(delay).await;
                }
            }
        }
    }
    
    Err(last_error.unwrap())
}

#[cfg(all(test, feature = "zarr"))]
mod tests {
    use super::*;
    use std::env;

    #[test]
    fn test_retry_config_delay_calculation() {
        let config = RetryConfig {
            max_retries: 3,
            base_delay_ms: 100,
            exponential_backoff: true,
            max_delay_ms: 30_000,
        };
        
        assert_eq!(config.calculate_delay(0), Duration::from_millis(0));
        assert_eq!(config.calculate_delay(1), Duration::from_millis(100));
        assert_eq!(config.calculate_delay(2), Duration::from_millis(200));
        assert_eq!(config.calculate_delay(3), Duration::from_millis(400));
    }
    
    #[test]
    fn test_retry_config_linear_backoff() {
        let config = RetryConfig {
            max_retries: 3,
            base_delay_ms: 100,
            exponential_backoff: false,
            max_delay_ms: 30_000,
        };
        
        assert_eq!(config.calculate_delay(0), Duration::from_millis(0));
        assert_eq!(config.calculate_delay(1), Duration::from_millis(100));
        assert_eq!(config.calculate_delay(2), Duration::from_millis(100));
        assert_eq!(config.calculate_delay(3), Duration::from_millis(100));
    }

    #[test]
    fn test_cloud_auth_env_var_loading() {
        // Set up test environment variables
        env::set_var("AWS_ACCESS_KEY_ID", "test_access_key");
        env::set_var("AWS_SECRET_ACCESS_KEY", "test_secret_key");
        env::set_var("AWS_REGION", "us-west-2");
        
        let provider = CloudAuthProvider::new(
            CloudAuth::default(),
            StreamingConfig::default(),
        );
        
        let auth = provider.load_auth().unwrap();
        
        assert_eq!(auth.aws_access_key_id, Some("test_access_key".to_string()));
        assert_eq!(auth.aws_secret_access_key, Some("test_secret_key".to_string()));
        assert_eq!(auth.aws_region, Some("us-west-2".to_string()));
        
        // Clean up
        env::remove_var("AWS_ACCESS_KEY_ID");
        env::remove_var("AWS_SECRET_ACCESS_KEY");
        env::remove_var("AWS_REGION");
    }

    #[test]
    fn test_config_override_env_vars() {
        // Set environment variable
        env::set_var("AWS_ACCESS_KEY_ID", "env_access_key");
        
        let mut auth_config = CloudAuth::default();
        auth_config.aws_access_key_id = Some("config_access_key".to_string());
        
        let provider = CloudAuthProvider::new(auth_config, StreamingConfig::default());
        let auth = provider.load_auth().unwrap();
        
        // Config should override env var
        assert_eq!(auth.aws_access_key_id, Some("config_access_key".to_string()));
        
        // Clean up
        env::remove_var("AWS_ACCESS_KEY_ID");
    }
}
