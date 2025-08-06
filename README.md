# QIBT Rust

Quasi-Isentropic Back-Trajectory Analysis Tool in Rust.

**Version 1.1.0** - Now with unified data format support, cloud integration, and advanced streaming capabilities!

## Build Instructions

Before building, ensure you have Rust and Cargo installed. You can install Rust using `rustup`:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

To build the project, run:

```bash
cargo build --release
```

This will produce the executable in the `target/release` directory.

## Run Instructions

To run the project, use the following command format:

```bash
./target/release/qibt_rust <SUBCOMMAND> [OPTIONS]
```

- **validate**: Validate Rust implementation against a Fortran reference.
- **benchmark**: Benchmark performance with different thread counts.
- **run**: Run trajectory computation.
- **test-sample**: Run Rust version on a 1-2 day sample for validation.

### Supported Formats:

The CLI now supports directories and URLs as inputs for Zarr and NetCDF formats. The input format can be specified or automatically detected:
- `--format netcdf` for NetCDF inputs.
- `--format zarr` for Zarr inputs.
- `--format auto` for automatic detection (default).

### Examples:

**Basic NetCDF file processing:**
```bash
./target/release/qibt_rust run --input data.nc --output trajectory.nc
```

**Processing a Zarr dataset:**
```bash
./target/release/qibt_rust run --input data.zarr --output trajectory.nc --format zarr
```

**Processing from a URL (with automatic format detection):**
```bash
./target/release/qibt_rust run --input https://example.com/data.zarr --output trajectory.nc
```

**Processing with custom parameters:**
```bash
./target/release/qibt_rust run \
    --input /path/to/data \
    --output trajectory.nc \
    --format auto \
    --parcels 500 \
    --threads 16 \
    --start-lat 40.0 \
    --start-lon -105.0 \
    --length 72
```

**Benchmarking with different formats:**
```bash
./target/release/qibt_rust benchmark \
    --input data.zarr \
    --format zarr \
    --parcels 1000 \
    --thread-counts 1,8,16,32
```

**Testing with format validation:**
```bash
./target/release/qibt_rust test-single-file \
    --file /path/to/dataset \
    --format auto \
    --variable temperature \
    --info
```

## Zarr Support

QIBT Rust includes optional Zarr support for reading chunked, compressed scientific datasets. To enable Zarr support:

### Build with Zarr:

```bash
cargo build --release --features zarr
```

### Zarr Features:

The following Zarr-related features are available:
- `zarr` - Basic Zarr support with local and HTTP access
- `zarr_s3` - Amazon S3 backend support
- `zarr_gcs` - Google Cloud Storage backend support
- `zarr_azure` - Azure Blob Storage backend support
- `zarr_async` - Asynchronous Zarr operations

### Zarr Input Examples:

**Local Zarr directory:**
```bash
./target/release/qibt_rust run --input /data/dataset.zarr --format zarr
```

**Zarr from HTTP URL:**
```bash
./target/release/qibt_rust run --input https://example.com/data.zarr --format zarr
```

**Zarr with automatic detection:**
```bash
./target/release/qibt_rust run --input ./weather_data.zarr --format auto
```

## ☁️ Cloud Storage Support

### Authentication Setup

QIBT Rust supports multiple cloud authentication methods:

#### AWS S3
```bash
# Method 1: Environment variables (recommended)
export AWS_ACCESS_KEY_ID=your_access_key
export AWS_SECRET_ACCESS_KEY=your_secret_key
export AWS_REGION=us-west-2

# Method 2: CLI parameters  
./target/release/qibt_rust run \
    --input s3://bucket/data.zarr \
    --aws-access-key-id your_key \
    --aws-secret-access-key your_secret \
    --aws-region us-west-2
```

#### Google Cloud Storage
```bash
# Service account key file
export GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account-key.json
export GOOGLE_CLOUD_PROJECT=your-project-id

# Use with GCS URLs
./target/release/qibt_rust run --input gs://bucket/data.zarr --output results.nc
```

#### Azure Blob Storage
```bash
# Storage account credentials
export AZURE_STORAGE_ACCOUNT=your_account
export AZURE_STORAGE_KEY=your_key

# Use with Azure URLs
./target/release/qibt_rust run \
    --input abfs://container/account.dfs.core.windows.net/data.zarr \
    --output results.nc
```

#### Custom S3-Compatible Services (MinIO, etc.)
```bash
export AWS_ENDPOINT_URL=https://minio.example.com
export AWS_ACCESS_KEY_ID=minio_access_key
export AWS_SECRET_ACCESS_KEY=minio_secret_key

./target/release/qibt_rust run --input s3://bucket/data.zarr --output results.nc
```

### Cloud Examples

**Processing large climate datasets from AWS:**
```bash
cargo build --release --features zarr_s3

./target/release/qibt_rust run \
    --input s3://climate-data-bucket/era5/temperature.zarr \
    --output local_results.nc \
    --chunk-size 16777216 \
    --max-concurrent-requests 8 \
    --parcels 1000
```

**Streaming analysis with optimization:**
```bash
./target/release/qibt_rust benchmark \
    --input https://data.example.com/weather.zarr \
    --format zarr \
    --thread-counts 4,8,16 \
    --chunk-size 8388608 \
    --request-timeout-secs 60
```

**Multi-cloud workflow:**
```bash
# Process data from Google Cloud, output to AWS
export GOOGLE_APPLICATION_CREDENTIALS=service-key.json
export AWS_ACCESS_KEY_ID=output_key
export AWS_SECRET_ACCESS_KEY=output_secret

./target/release/qibt_rust run \
    --input gs://weather-source/model-data.zarr \
    --output s3://results-bucket/trajectory.nc
```

## Dependencies

This project uses several dependencies. They are automatically installed when using Cargo:

- `netcdf`
- `ndarray`
- `rayon`
- `chrono`
- `clap`
- `num-traits`
- `thiserror`
- `rand`
- `crossbeam-channel`
- `flamegraph` (optional)
- `zarrs` (optional, with zarr feature)
- `zarrs_object_store` (optional, with zarr feature)
- `tokio` (optional, with zarr feature)

Intel MKL is not required.

## Example Job Script

Below is an example shell script to run the compiled binary:

```bash
#!/bin/bash

# Example job script for QIBT Rust

./target/release/qibt_rust run \
    --input test_data/wrfout_d01_2023-07-31_12:00:00.nc \
    --output results/output_trajectory.nc \
    --parcels 100 \
    --threads 8
```

Make sure the input data file exists at the specified path.

## 🚀 Performance Optimization

### Memory and Network Tuning

**For large local datasets:**
```bash
# Increase memory usage for better performance
./target/release/qibt_rust run \
    --input large_dataset.zarr \
    --chunk-size 33554432 \    # 32MB chunks
    --threads 16
```

**For cloud datasets with limited bandwidth:**
```bash
# Optimize for slower connections
./target/release/qibt_rust run \
    --input s3://bucket/data.zarr \
    --chunk-size 4194304 \     # 4MB chunks
    --max-concurrent-requests 2 \
    --request-timeout-secs 120
```

**For high-bandwidth cloud access:**
```bash
# Maximize throughput
./target/release/qibt_rust run \
    --input gs://bucket/data.zarr \
    --chunk-size 67108864 \    # 64MB chunks  
    --max-concurrent-requests 16 \
    --request-timeout-secs 30
```

### Benchmarking Your Setup

```bash
# Test different configurations
./target/release/qibt_rust benchmark \
    --input your_dataset.zarr \
    --thread-counts 1,4,8,16,32 \
    --chunk-size 8388608 \
    --parcels 500

# Compare formats
./target/release/qibt_rust benchmark --input data.nc --thread-counts 4,8,16
./target/release/qibt_rust benchmark --input data.zarr --thread-counts 4,8,16 --format zarr
```

## 🔧 Advanced Usage

### Feature Combinations

```bash
# Basic Zarr support
cargo build --release --features zarr

# All cloud backends
cargo build --release --features zarr,zarr_s3,zarr_gcs,zarr_azure

# With async support for large-scale processing
cargo build --release --features zarr,zarr_s3,zarr_async

# Development build with profiling
cargo build --features zarr,optimized
```

### Custom Configuration Examples

**Research cluster with shared storage:**
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128GB

./target/release/qibt_rust run \
    --input /shared/climate_data/era5.zarr \
    --output $SCRATCH/results/trajectory_$(date +%Y%m%d).nc \
    --threads 32 \
    --parcels 10000 \
    --chunk-size 134217728   # 128MB for high-memory system
```

## ⚠️ Troubleshooting

### Common Issues

**"Zarr support not compiled in"**
```bash
# Solution: Rebuild with Zarr features
cargo clean
cargo build --release --features zarr
```

**"Unable to detect format"**
```bash
# Solution 1: Specify format explicitly
./target/release/qibt_rust run --input data.zarr --format zarr

# Solution 2: Check file structure
ls -la data.zarr/  # Should contain .zgroup or .zarray files

# Solution 3: Verify file accessibility
file data.nc  # Should show NetCDF format
```

**Performance issues with cloud data**
```bash
# Reduce concurrent requests
./target/release/qibt_rust run \
    --input s3://bucket/data.zarr \
    --max-concurrent-requests 2 \
    --request-timeout-secs 120

# Increase chunk size for better throughput
./target/release/qibt_rust run \
    --input gs://bucket/data.zarr \
    --chunk-size 33554432  # 32MB
```

**Memory issues with large datasets**
```bash
# Reduce chunk size to lower memory usage
./target/release/qibt_rust run \
    --input large_dataset.zarr \
    --chunk-size 2097152 \     # 2MB chunks
    --threads 4                # Fewer threads

# Use streaming mode for very large datasets
./target/release/qibt_rust run \
    --input https://data.example.com/huge.zarr \
    --chunk-size 1048576 \     # 1MB chunks
    --max-concurrent-requests 1
```

### Debug Information

```bash
# Enable verbose logging
RUST_LOG=debug ./target/release/qibt_rust run --input data.zarr

# Get detailed error information
RUST_BACKTRACE=1 ./target/release/qibt_rust run --input problematic_data.zarr

# Test dataset accessibility
./target/release/qibt_rust test-single-file \
    --file data.zarr \
    --format auto \
    --info
```

## 📚 Examples and Tutorials

See the `examples/` directory for comprehensive examples:

- `cloud_streaming_demo.rs` - Cloud authentication and streaming patterns
- `unified_loader_demo.rs` - Format detection and unified interface usage  
- `unified_trajectory_demo.rs` - Trajectory computation with different readers

Run examples with:
```bash
# Cloud streaming demo
cargo run --example cloud_streaming_demo --features zarr_s3

# Unified loader demo
cargo run --example unified_loader_demo --features zarr

# Trajectory computation demo
cargo run --example unified_trajectory_demo --features zarr
```

## 🌐 Supported Data Sources

| Source Type | URL Format | Authentication | Build Features |
|-------------|------------|----------------|----------------|
| Local NetCDF | `/path/to/data.nc` | None | Default |
| Local Zarr | `/path/to/data.zarr` | None | `zarr` |
| HTTP/HTTPS | `https://example.com/data.zarr` | None | `zarr` |
| AWS S3 | `s3://bucket/data.zarr` | AWS credentials | `zarr_s3` |
| Google Cloud | `gs://bucket/data.zarr` | GCP service account | `zarr_gcs` |
| Azure Blob | `abfs://container/account.dfs.core.windows.net/data.zarr` | Azure credentials | `zarr_azure` |
| MinIO/Custom | `s3://bucket/data.zarr` + endpoint | S3-compatible credentials | `zarr_s3` |

## 🚀 What's Next

Upcoming features in future releases:
- Enhanced parallel streaming
- Built-in data validation tools
- Interactive trajectory visualization
- GPU support