# Migration Guide: QIBT Rust 1.0.0 → 1.1.0

This guide helps you migrate from QIBT Rust 1.0.0 to 1.1.0, which introduces significant new capabilities while maintaining full backward compatibility.

## Overview

Version 1.1.0 introduces:
- Unified DataReader interface supporting multiple formats
- Zarr format support with cloud integration
- Advanced streaming capabilities
- Enhanced CLI options

**Good news**: All existing 1.0.0 code continues to work without changes!

## What's New

### 1. Unified Data Format Support

#### Before (1.0.0) - NetCDF Only
```bash
# Only supported NetCDF files
./qibt_rust run --input data.nc --output trajectory.nc
```

#### After (1.1.0) - Multiple Formats
```bash
# NetCDF (unchanged)
./qibt_rust run --input data.nc --output trajectory.nc

# Zarr datasets (new!)
./qibt_rust run --input data.zarr --output trajectory.nc --format zarr

# Automatic format detection (new!)
./qibt_rust run --input data.zarr --output trajectory.nc --format auto

# Cloud datasets (new!)
./qibt_rust run --input s3://bucket/data.zarr --output trajectory.nc
```

### 2. Enhanced Build Options

#### Before (1.0.0)
```bash
# Basic build
cargo build --release
```

#### After (1.1.0) - Feature Flags
```bash
# Basic build (NetCDF only, same as before)
cargo build --release

# With Zarr support
cargo build --release --features zarr

# With cloud backends
cargo build --release --features zarr_s3
cargo build --release --features zarr_gcs  
cargo build --release --features zarr_azure

# All features
cargo build --release --features zarr,zarr_s3,zarr_gcs,zarr_azure,zarr_async
```

### 3. New CLI Options

All existing CLI options remain the same. New options include:

```bash
# Format specification
--format auto|netcdf|zarr    # Default: auto

# Cloud authentication (when using cloud URLs)
--aws-access-key-id KEY
--aws-secret-access-key SECRET  
--aws-region REGION
--gcp-service-account-key PATH
--azure-storage-account ACCOUNT

# Streaming configuration
--chunk-size BYTES           # Default: 8388608 (8MB)
--max-concurrent-requests N  # Default: 4
--request-timeout-secs N     # Default: 30
```

## Migration Scenarios

### Scenario 1: Existing NetCDF Workflows

**No changes required!** Your existing scripts continue to work:

```bash
# This continues to work exactly as before
./qibt_rust run --input weather.nc --output results.nc --parcels 100 --threads 8
./qibt_rust benchmark --input data.nc --thread-counts 1,4,8,16
./qibt_rust validate --input reference.nc
```

### Scenario 2: Adding Zarr Support to Existing Workflows

#### Step 1: Rebuild with Zarr Support
```bash
cargo build --release --features zarr
```

#### Step 2: Use Same Commands with Zarr Data
```bash
# Replace .nc with .zarr in your existing commands
./qibt_rust run --input weather.zarr --output results.nc --parcels 100 --threads 8
./qibt_rust benchmark --input data.zarr --thread-counts 1,4,8,16 --format zarr
```

### Scenario 3: Migrating to Cloud Data Sources

#### Step 1: Set Up Authentication
```bash
# Option A: Environment variables (recommended)
export AWS_ACCESS_KEY_ID=your_access_key
export AWS_SECRET_ACCESS_KEY=your_secret_key
export AWS_REGION=us-west-2

# Option B: Use CLI parameters
# (credentials passed via command line)
```

#### Step 2: Build with Cloud Support
```bash
cargo build --release --features zarr_s3
```

#### Step 3: Update Data Sources
```bash
# Before: local file
./qibt_rust run --input /local/path/data.nc --output results.nc

# After: cloud dataset
./qibt_rust run --input s3://my-bucket/data.zarr --output results.nc
```

### Scenario 4: Optimizing for Large Datasets

Take advantage of new streaming capabilities:

```bash
# Optimize chunk size for your network/memory constraints
./qibt_rust run \
    --input s3://weather-data/large-dataset.zarr \
    --output results.nc \
    --chunk-size 16777216 \     # 16MB chunks
    --max-concurrent-requests 8 \
    --request-timeout-secs 60
```

## Code Migration (for library users)

If you're using QIBT Rust as a library, the DataReader trait provides new capabilities:

### Before (1.0.0) - Direct NetCDF Usage
```rust
use qibt_rust::data_io::NetCDFReader;

let mut reader = NetCDFReader::new("data.nc");
// NetCDF-specific code...
```

### After (1.1.0) - Unified Interface (Optional)
```rust
use qibt_rust::io::{Dataset, DataReader};

// Option 1: Use the unified Dataset interface (recommended)
let dataset = Dataset::open("data.nc")?;  // Works with NetCDF or Zarr
let variables = dataset.list_variables()?;

// Option 2: Use specific readers (when you need format-specific features)
let mut netcdf_reader = NetCDFReader::new("data.nc");
let zarr_reader = create_zarr_reader("data.zarr")?;

// Option 3: Generic code that works with any format
fn process_dataset<T: DataReader>(reader: &T) -> Result<(), Error> {
    let vars = reader.list_variables()?;
    // ... same code works for NetCDF, Zarr, or future formats
    Ok(())
}
```

## Performance Considerations

### Memory Usage
- **1.0.0**: All data loaded into memory at once
- **1.1.0**: Streaming support reduces memory usage for large datasets
  - Configure `--chunk-size` based on available memory
  - Default 8MB chunks work well for most systems

### Network Usage  
- **New in 1.1.0**: Range requests minimize data transfer
- Only requested data portions are downloaded
- Intelligent caching reduces repeated requests

### Parallel Processing
- **1.0.0**: File-level parallelism
- **1.1.0**: Adds chunk-level parallelism for cloud data
  - Configure `--max-concurrent-requests` based on network capacity
  - Higher values = faster download but more bandwidth usage

## Environment Variables

### New in 1.1.0
```bash
# AWS Authentication
AWS_ACCESS_KEY_ID / QIBT_AWS_ACCESS_KEY_ID
AWS_SECRET_ACCESS_KEY / QIBT_AWS_SECRET_ACCESS_KEY  
AWS_SESSION_TOKEN / QIBT_AWS_SESSION_TOKEN
AWS_REGION / QIBT_AWS_REGION
AWS_ENDPOINT_URL / QIBT_ENDPOINT_URL

# Google Cloud Authentication
GOOGLE_APPLICATION_CREDENTIALS / QIBT_GCP_SERVICE_ACCOUNT_KEY
GOOGLE_CLOUD_PROJECT / QIBT_GCP_PROJECT_ID

# Azure Authentication  
AZURE_STORAGE_ACCOUNT / QIBT_AZURE_STORAGE_ACCOUNT
AZURE_STORAGE_KEY / QIBT_AZURE_STORAGE_KEY
```

### Existing (Still Supported)
All existing environment variables continue to work as before.

## Troubleshooting

### "Zarr support not compiled in"
**Solution**: Rebuild with Zarr features:
```bash
cargo build --release --features zarr
```

### "Unable to detect format"
**Solutions**:
1. Specify format explicitly: `--format zarr` or `--format netcdf`
2. Check file exists and is readable
3. For Zarr directories, ensure `.zgroup` or `.zarray` files are present

### Cloud Authentication Errors
**Solutions**:
1. Verify environment variables are set correctly
2. Check credential permissions for the specific bucket/container
3. Use CLI parameters to override environment variables
4. For custom endpoints, set `QIBT_ENDPOINT_URL`

### Performance Issues with Cloud Data
**Solutions**:
1. Increase chunk size: `--chunk-size 16777216` (16MB)
2. Reduce concurrent requests: `--max-concurrent-requests 2`
3. Increase timeout: `--request-timeout-secs 120`
4. Check network bandwidth and latency

## Testing Your Migration

### 1. Verify Backward Compatibility
```bash
# These should work exactly as before
./qibt_rust run --input old_data.nc --output test.nc
./qibt_rust benchmark --input old_data.nc --thread-counts 1,4
```

### 2. Test New Features
```bash
# Test format detection
./qibt_rust run --input data.zarr --output test.nc --format auto

# Test cloud access (if applicable)
./qibt_rust run --input s3://test-bucket/data.zarr --output test.nc
```

### 3. Performance Comparison
```bash
# Compare performance between versions
./qibt_rust benchmark --input large_dataset.nc --thread-counts 1,4,8,16
./qibt_rust benchmark --input large_dataset.zarr --thread-counts 1,4,8,16 --features zarr
```

## Getting Help

- **Documentation**: All existing documentation applies to NetCDF workflows
- **Examples**: Check `examples/` directory for new usage patterns
- **Issues**: Report problems with migration on GitHub
- **Performance**: Use `benchmark` subcommand to optimize settings

## Summary

The migration to 1.1.0 is **opt-in**:
- Existing workflows continue unchanged
- New features available when explicitly enabled
- No breaking changes to existing APIs
- Gradual adoption supported

This approach ensures you can upgrade safely while taking advantage of new capabilities at your own pace.
