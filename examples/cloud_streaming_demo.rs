/// Example demonstrating cloud authentication and streaming capabilities
/// 
/// This example shows how to:
/// 1. Configure cloud authentication from environment variables or CLI
/// 2. Stream large datasets from cloud storage without full download
/// 3. Use ranged GETs for efficient data access
/// 4. Implement retry logic with exponential backoff
/// 
/// Run with: cargo run --example cloud_streaming_demo --features zarr_s3

use qibt_rust::config::{CloudAuth, StreamingConfig};
use qibt_rust::io::cloud_auth::{CloudAuthProvider, CloudAuthError};
use qibt_rust::io::streaming_reader::{StreamingReader, DataRange, StreamingError};

#[cfg(feature = "zarr")]
use std::env;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== QIBT Cloud Streaming Demo ===\n");

    #[cfg(not(feature = "zarr"))]
    {
        println!("This example requires the 'zarr' feature to be enabled.");
        println!("Run with: cargo run --example cloud_streaming_demo --features zarr_s3");
        return Ok(());
    }

    #[cfg(feature = "zarr")]
    {
        demo_environment_variables()?;
        demo_configuration_structures().await?;
        demo_streaming_concepts().await?;
        demo_retry_logic().await?;
        demo_cloud_urls().await?;
    }

    Ok(())
}

#[cfg(feature = "zarr")]
fn demo_environment_variables() -> Result<(), CloudAuthError> {
    println!("1. Environment Variable Configuration");
    println!("=====================================");
    
    // Show supported environment variables
    println!("Supported environment variables:");
    println!("  AWS:");
    println!("    AWS_ACCESS_KEY_ID / QIBT_AWS_ACCESS_KEY_ID");
    println!("    AWS_SECRET_ACCESS_KEY / QIBT_AWS_SECRET_ACCESS_KEY");
    println!("    AWS_SESSION_TOKEN / QIBT_AWS_SESSION_TOKEN");
    println!("    AWS_REGION / QIBT_AWS_REGION");
    println!("  GCP:");
    println!("    GOOGLE_APPLICATION_CREDENTIALS / QIBT_GCP_SERVICE_ACCOUNT_KEY");
    println!("    GOOGLE_CLOUD_PROJECT / QIBT_GCP_PROJECT_ID");
    println!("  Azure:");
    println!("    AZURE_STORAGE_ACCOUNT / QIBT_AZURE_STORAGE_ACCOUNT");
    println!("    AZURE_STORAGE_KEY / QIBT_AZURE_STORAGE_KEY");
    println!("  Custom:");
    println!("    QIBT_ENDPOINT_URL / AWS_ENDPOINT_URL\n");

    // Create auth provider and load from environment
    let auth_provider = CloudAuthProvider::new(
        CloudAuth::default(),
        StreamingConfig::default(),
    );
    
    let auth = auth_provider.load_auth()?;
    
    println!("Current environment configuration:");
    println!("  AWS Access Key ID: {}", 
        auth.aws_access_key_id.as_deref().unwrap_or("Not set"));
    println!("  AWS Region: {}", 
        auth.aws_region.as_deref().unwrap_or("Not set"));
    println!("  GCP Service Account: {}", 
        auth.gcp_service_account_key.as_deref().unwrap_or("Not set"));
    println!("  Custom Endpoint: {}", 
        auth.endpoint_url.as_deref().unwrap_or("Not set"));
    
    println!();
    Ok(())
}

#[cfg(feature = "zarr")]
async fn demo_configuration_structures() -> Result<(), StreamingError> {
    println!("2. Configuration Structures");
    println!("===========================");
    
    // Demonstrate cloud authentication configuration
    let mut cloud_auth = CloudAuth::default();
    cloud_auth.aws_access_key_id = Some("demo-access-key".to_string());
    cloud_auth.aws_region = Some("us-west-2".to_string());
    cloud_auth.endpoint_url = Some("https://minio.example.com".to_string());
    
    println!("CloudAuth configuration:");
    println!("  {:?}", cloud_auth);
    
    // Demonstrate streaming configuration
    let streaming_config = StreamingConfig {
        chunk_size: 8 * 1024 * 1024,     // 8MB chunks
        max_concurrent_requests: 5,       // Conservative for demo
        request_timeout_secs: 60,         // 1 minute timeout
        max_retries: 5,                   // More retries for demo
        retry_delay_ms: 200,              // 200ms base delay
        exponential_backoff: true,        // Use exponential backoff
        buffer_size: 32 * 1024 * 1024,   // 32MB buffer
        enable_range_requests: true,      // Enable range requests
    };
    
    println!("\nStreamingConfig:");
    println!("  Chunk size: {} MB", streaming_config.chunk_size / (1024 * 1024));
    println!("  Max concurrent requests: {}", streaming_config.max_concurrent_requests);
    println!("  Request timeout: {}s", streaming_config.request_timeout_secs);
    println!("  Max retries: {}", streaming_config.max_retries);
    println!("  Exponential backoff: {}", streaming_config.exponential_backoff);
    println!("  Range requests: {}", streaming_config.enable_range_requests);
    
    println!();
    Ok(())
}

#[cfg(feature = "zarr")]
async fn demo_streaming_concepts() -> Result<(), StreamingError> {
    println!("3. Streaming Concepts");
    println!("====================");
    
    // Create a streaming reader
    let auth_provider = CloudAuthProvider::new(
        CloudAuth::default(),
        StreamingConfig::default(),
    );
    
    let reader = StreamingReader::new(
        StreamingConfig::default(),
        auth_provider,
    );
    
    // Demonstrate DataRange usage
    println!("DataRange examples:");
    
    let range1 = DataRange::new(0, Some(1024));
    println!("  Fixed range: bytes {}-{} ({})", 
        range1.start, 
        range1.end.unwrap_or(0), 
        range1.length().unwrap_or(0));
    
    let range2 = DataRange::from_offset_and_length(2048, 512);
    println!("  Offset+length: bytes {}-{} ({} bytes)", 
        range2.start, 
        range2.end.unwrap_or(0), 
        range2.length().unwrap_or(0));
    
    let range3 = DataRange::new(1024, None);
    println!("  Open-ended: bytes {}- (to end)", range3.start);
    
    // Show cache statistics
    let (entries, size) = reader.get_cache_stats().await;
    println!("\nCache statistics:");
    println!("  Cached entries: {}", entries);
    println!("  Cache size: {} bytes", size);
    
    // Demonstrate multiple concurrent ranges
    let ranges = vec![
        DataRange::new(0, Some(1024)),
        DataRange::new(1024, Some(2048)),
        DataRange::new(2048, Some(3072)),
    ];
    
    println!("\nConcurrent range requests ({} ranges):", ranges.len());
    for (i, range) in ranges.iter().enumerate() {
        println!("  Range {}: {}-{}", i, range.start, range.end.unwrap_or(0));
    }
    
    println!();
    Ok(())
}

#[cfg(feature = "zarr")]
async fn demo_retry_logic() -> Result<(), StreamingError> {
    println!("4. Retry Logic with Exponential Backoff");
    println!("=======================================");
    
    use qibt_rust::io::cloud_auth::RetryConfig;
    use std::time::Duration;
    
    // Demonstrate retry configuration
    let streaming_config = StreamingConfig {
        max_retries: 4,
        retry_delay_ms: 100,
        exponential_backoff: true,
        ..Default::default()
    };
    
    let retry_config = RetryConfig::from(&streaming_config);
    
    println!("Retry configuration:");
    println!("  Max retries: {}", retry_config.max_retries);
    println!("  Base delay: {}ms", retry_config.base_delay_ms);
    println!("  Exponential backoff: {}", retry_config.exponential_backoff);
    println!("  Max delay: {}ms", retry_config.max_delay_ms);
    
    println!("\nRetry delay progression:");
    for attempt in 0..=retry_config.max_retries {
        let delay = retry_config.calculate_delay(attempt);
        println!("  Attempt {}: {}ms delay", attempt, delay.as_millis());
    }
    
    // Demonstrate linear backoff
    let linear_config = RetryConfig {
        max_retries: 4,
        base_delay_ms: 200,
        exponential_backoff: false,
        max_delay_ms: 30_000,
    };
    
    println!("\nLinear backoff progression:");
    for attempt in 0..=linear_config.max_retries {
        let delay = linear_config.calculate_delay(attempt);
        println!("  Attempt {}: {}ms delay", attempt, delay.as_millis());
    }
    
    println!();
    Ok(())
}

#[cfg(feature = "zarr")]
async fn demo_cloud_urls() -> Result<(), StreamingError> {
    println!("5. Cloud URL Examples");
    println!("====================");
    
    // Show examples of supported URL schemes
    let example_urls = vec![
        ("AWS S3", "s3://my-bucket/path/to/data.zarr"),
        ("Google Cloud Storage", "gs://my-bucket/path/to/data.zarr"),
        ("Azure Blob Storage", "abfs://container/account.dfs.core.windows.net/path/data.zarr"),
        ("HTTP/HTTPS", "https://example.com/data.zarr"),
        ("Custom S3 endpoint", "s3://bucket/data.zarr (with QIBT_ENDPOINT_URL=https://minio.example.com)"),
    ];
    
    println!("Supported URL schemes:");
    for (provider, url) in example_urls {
        println!("  {}: {}", provider, url);
    }
    
    // Demonstrate URL parsing logic
    println!("\nURL scheme detection:");
    let test_urls = vec![
        "s3://weather-data/zarr/temperature.zarr",
        "gs://climate-bucket/humidity.zarr", 
        "https://data.example.com/wind.zarr",
        "abfs://container/account.dfs.core.windows.net/pressure.zarr",
    ];
    
    for url in test_urls {
        let scheme = if url.starts_with("s3://") {
            "AWS S3"
        } else if url.starts_with("gs://") {
            "Google Cloud Storage"
        } else if url.starts_with("abfs://") {
            "Azure Blob Storage"
        } else if url.starts_with("http") {
            "HTTP/HTTPS"
        } else {
            "Unknown"
        };
        println!("  {}: {} scheme", url, scheme);
    }
    
    println!();
    Ok(())
}

#[cfg(feature = "zarr")]
async fn demo_usage_patterns() -> Result<(), StreamingError> {
    println!("6. Usage Patterns");
    println!("================");
    
    println!("Example usage patterns:");
    
    println!("\n// 1. Stream processing with chunked reads");
    println!("let reader = StreamingReader::new(config, auth_provider);");
    println!("reader.initialize(\"s3://bucket/data.zarr\").await?;");
    println!("reader.read_chunked(\"temperature\", |chunk, offset| {{");
    println!("    // Process chunk without loading entire dataset");
    println!("    process_temperature_chunk(chunk, offset)");
    println!("}}).await?;");
    
    println!("\n// 2. Selective data access with range requests");
    println!("let ranges = vec![");
    println!("    DataRange::new(0, Some(1024*1024)),      // First MB");
    println!("    DataRange::new(1024*1024*100, None),     // Skip to 100MB");
    println!("];");
    println!("let data = reader.read_multiple_ranges(\"humidity\", ranges).await?;");
    
    println!("\n// 3. Metadata inspection without data download");
    println!("let metadata = reader.get_object_metadata(\"wind.zarr\").await?;");
    println!("println!(\"Size: {} bytes, Modified: {{}}\", metadata.size, metadata.last_modified);");
    
    println!("\n// 4. Environment-based authentication");
    println!("export AWS_ACCESS_KEY_ID=your_key");
    println!("export AWS_SECRET_ACCESS_KEY=your_secret");
    println!("export AWS_REGION=us-west-2");
    println!("# Authentication is automatic from environment");
    
    println!("\n// 5. CLI override of environment variables");
    println!("./qibt_rust benchmark --input s3://bucket/data.zarr \\");
    println!("  --aws-access-key-id override_key \\");
    println!("  --aws-region us-east-1 \\");
    println!("  --chunk-size 8388608 \\");
    println!("  --max-concurrent-requests 5");
    
    println!();
    Ok(())
}
