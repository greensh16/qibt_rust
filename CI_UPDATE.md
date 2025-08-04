# CI and Feature Gates Update

## Summary

This update implements the following CI and feature gate improvements:

### 1. ✅ Zarr feature enabled in default test matrix

- Updated `.github/workflows/ci.yml` to include `zarr` and `zarr,optimized` in the default test matrix
- CI now tests zarr functionality by default on all platforms (Ubuntu, Windows, macOS)
- Both stable and beta Rust versions are tested with zarr features

### 2. ✅ Added job caching for large sample data

- Implemented caching strategy for test data in `~/.cache/qibt_test_data/`
- Cache uses file checksums for invalidation (`test_data/**/*.md5`, `test_data/**/*.sha256`)
- Downloads are only performed if cached data is not available
- Significant time savings for CI runs with large datasets

### 3. ✅ Verified NetCDF-only builds work without zarr

- Added specific test jobs for `--no-default-features --features netcdf`
- Implemented proper conditional compilation with `#[cfg(feature = "zarr")]`
- Created dedicated CI job `feature-compatibility` that verifies:
  - Minimal NetCDF build works without zarr dependencies
  - Individual zarr features can be enabled
  - Various feature combinations work correctly
  - Dependency tree doesn't include zarr when not enabled

## Feature Organization

### Available Features

```toml
[features]
default = []
netcdf = []                    # Basic NetCDF support
optimized = ["flamegraph"]     # Performance profiling
f32_precision = []             # Single precision floats

# Zarr support with cloud backends
zarr = [                       # Core zarr functionality
    "dep:zarrs",
    "dep:zarrs_object_store", 
    "dep:zarrs_http",
    "dep:object_store",
    "dep:async-compression",
    "dep:tokio",
    "dep:serde_json",
    "dep:bytes",
    "dep:futures",
    "object_store/http",
    "object_store/fs"
]

# Individual cloud backend features (require zarr feature)
zarr_s3 = ["zarr", "zarrs_object_store/aws", "object_store/aws"]
zarr_gcs = ["zarr", "zarrs_object_store/gcp", "object_store/gcp"] 
zarr_http = ["zarr", "zarrs_object_store/http", "object_store/http"]
zarr_azure = ["zarr", "zarrs_object_store/azure", "object_store/azure"]
zarr_async = ["zarr", "zarrs/async"]
```

### Build Examples

```bash
# NetCDF-only build (minimal dependencies)
cargo build --no-default-features --features netcdf

# Default zarr build
cargo build --features zarr

# Zarr with cloud backends
cargo build --features "zarr,zarr_s3,zarr_gcs"

# Optimized build with zarr
cargo build --features "zarr,optimized"
```

## CI Workflow Structure

### Test Matrix
- **Operating Systems**: Ubuntu, Windows, macOS
- **Rust Versions**: stable, beta
- **Default Features**: `zarr`, `zarr,optimized`
- **NetCDF-only**: Special test cases for `--no-default-features --features netcdf`

### Jobs

1. **Test Suite** - Main testing with zarr features enabled
2. **Clippy** - Linting for both zarr and netcdf-only builds  
3. **Format** - Code formatting checks
4. **Documentation** - Doc generation for all feature combinations
5. **Feature Compatibility** - Comprehensive feature combination testing

### Data Caching

The CI implements a two-tier caching strategy:

1. **Rust Dependencies**: Standard cargo cache for faster compilation
2. **Large Sample Data**: Custom cache for test datasets with intelligent invalidation

```yaml
- name: Cache large sample data
  uses: actions/cache@v4
  with:
    path: |
      test_data/
      ~/.cache/qibt_test_data/
    key: ${{ runner.os }}-test-data-v2-${{ hashFiles('test_data/**/*.md5', 'test_data/**/*.sha256') }}
```

## Verification

Run the verification script to test locally:

```bash
./verify_ci.sh
```

This script validates:
- NetCDF-only builds work without zarr dependencies
- Feature gates are properly implemented
- No zarr dependencies leak into netcdf-only builds
- All builds complete successfully

## Migration Notes

### For Users Who Don't Need Zarr

The default build behavior remains unchanged. Users who only need NetCDF functionality can build with:

```bash
cargo build --no-default-features --features netcdf
```

This ensures minimal dependencies and faster compilation times.

### For Zarr Users

Zarr functionality is now properly feature-gated and included in CI testing. Enable with:

```bash
cargo build --features zarr
```

Additional cloud backends can be enabled as needed:

```bash
cargo build --features "zarr,zarr_s3"  # AWS S3 support
cargo build --features "zarr,zarr_gcs" # Google Cloud Storage
```

## Files Modified

- `.github/workflows/ci.yml` - New comprehensive CI workflow
- `Cargo.toml` - Added `netcdf` feature, improved zarr feature organization
- `src/config.rs` - Added conditional compilation for zarr-specific fields
- `src/io/mod.rs` - Made zarr modules conditional
- `src/data_io/reader.rs` - Fixed naming conflicts
- `verify_ci.sh` - Local verification script
- `CI_UPDATE.md` - This documentation
