# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2024-12-28

### Added

#### Major Features
- **Unified DataReader Interface**: Generic trait-based system supporting multiple data formats
  - Common interface for NetCDF and Zarr datasets
  - Automatic format detection based on file structure and magic bytes
  - Extensible design for future format support (HDF5, GRIB, etc.)
  
- **Zarr Support with Cloud Integration**: 
  - Full Zarr v2 specification support with chunked, compressed arrays
  - Multi-cloud backend support (AWS S3, Google Cloud Storage, Azure Blob Storage)
  - HTTP/HTTPS streaming for remote datasets
  - Configurable cloud authentication with environment variable support
  - Streaming reader with intelligent caching and range requests
  
- **Advanced Cloud Authentication**:
  - Environment-based authentication (AWS_*, GOOGLE_*, AZURE_* variables)
  - Custom endpoint support for MinIO and S3-compatible services
  - CLI override capabilities for deployment flexibility
  - Secure credential handling with no plaintext exposure

#### Performance & Scalability
- **Intelligent Streaming**: 
  - Chunked data access with configurable buffer sizes
  - Range request optimization for selective data loading
  - Exponential backoff retry logic with configurable timeouts
  - Multi-threaded concurrent requests with rate limiting
  
- **Memory-Efficient Processing**:
  - Stream processing without full dataset download
  - Configurable chunk sizes (default 8MB) for memory optimization
  - Smart caching with LRU eviction policies
  - Zero-copy operations where possible

#### Mathematical & Physics Features
- **Advanced Interpolation**: Generic trait-based field data access
- **Physics Calculations**: Comprehensive atmospheric trajectory physics
- **Time Series Support**: DateTime-aware data readers with closest-time algorithms
- **Parallel Processing**: Multi-threaded trajectory computation with work stealing

#### Developer Experience
- **Comprehensive Examples**: 
  - `cloud_streaming_demo.rs` - Cloud authentication and streaming patterns
  - `unified_loader_demo.rs` - Format detection and unified interface usage
  - `unified_trajectory_demo.rs` - Trajectory computation with different readers
  
- **Feature Flags**: Modular compilation with optional Zarr support
  - `zarr` - Basic Zarr support with local and HTTP access
  - `zarr_s3` - Amazon S3 backend
  - `zarr_gcs` - Google Cloud Storage backend  
  - `zarr_azure` - Azure Blob Storage backend
  - `zarr_async` - Asynchronous operations
  
- **CLI Enhancements**:
  - `--format` parameter with auto-detection (`auto`, `netcdf`, `zarr`)
  - Cloud authentication parameters (`--aws-access-key-id`, etc.)
  - Streaming configuration (`--chunk-size`, `--max-concurrent-requests`)

### Technical Improvements
- **Error Handling**: Comprehensive error types with detailed context
- **Type Safety**: Generic trait system eliminates format-specific code duplication
- **Memory Safety**: Rust's ownership model prevents data races and memory leaks
- **Async Support**: Full async/await support for cloud operations
- **Configuration**: Structured configuration with builder patterns

### Dependencies
- Added `zarrs` (0.21.2) for Zarr format support
- Added `object_store` (0.12.3) for cloud backend abstraction
- Added `tokio` (1.47) for async runtime
- Added `futures` (0.3) for async stream processing
- Added `bytes` (1.5) for efficient byte buffer handling

### Breaking Changes
- None - All changes are backward compatible additions

### Migration Guide
- Existing NetCDF-only workflows continue to work unchanged
- To enable Zarr support: `cargo build --features zarr`
- For cloud support: `cargo build --features zarr_s3`
- No API changes required for existing code

## [1.0.0] - 2024-12-27

### Added
- Initial release of QIBT Rust
- Basic NetCDF trajectory computation
- Parallel processing capabilities
- Command-line interface with basic options
- Mathematical interpolation and physics modules
- Time utilities and series processing
- Benchmark framework for performance testing

### Features
- Quasi-Isentropic Back-Trajectory Analysis
- Multi-threaded trajectory computation
- NetCDF input/output support
- Comprehensive physics calculations for atmospheric dynamics
- Flexible configuration system
- Performance benchmarking tools

### Dependencies
- `netcdf` (0.11) for NetCDF file handling
- `ndarray` (0.15) with BLAS and Rayon support
- `rayon` (1.7) for parallel processing
- `chrono` (0.4) for time handling
- `clap` (4) for command-line parsing
- Supporting utilities: `num-traits`, `thiserror`, `rand`, `crossbeam-channel`

[1.1.0]: https://github.com/user/qibt_rust/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/user/qibt_rust/releases/tag/v1.0.0
