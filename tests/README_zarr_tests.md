# Zarr Reader Tests

This directory contains comprehensive unit and integration tests for the ZarrReader functionality.

## Test Organization

### Unit Tests (`unit_tests` module)
- Basic ZarrReader creation and configuration
- Variable name mapping functionality  
- Opening/closing operations
- Error handling for unsupported operations (when zarr feature disabled)

### Integration Tests (`integration_tests` module)
- Reading 4D variables from synthetic Zarr stores
- Variable slicing with stride support
- Multi-variable reading
- Coordinate variable reading
- Attribute access (global and variable-level)
- Value equality verification against known arrays

### Parity Tests (`parity_tests` module)
- Comparing NetCDF and Zarr outputs for identical datasets
- Metadata comparison between formats
- Variable info consistency
- Data value verification

### CLI Smoke Tests (`cli_smoke_tests` module)
- Remote ERA5 Zarr data access (requires `--ignored` flag)
- Auto-format detection testing
- CLI Zarr feature verification

### Mock Cloud Tests (`mock_cloud_tests` module)
- Mock S3, GCS, and HTTP backends for fast CI
- Cloud authentication configuration testing
- Simulated cloud object store operations

### Performance Tests (`performance_tests` module)
- Large dataset handling
- Concurrent access patterns
- Metadata access performance benchmarks

## Running the Tests

### Prerequisites

1. **Install test dependencies:**
   ```bash
   # For Python-based data generation
   pip install numpy zarr netcdf4
   ```

2. **Generate test data (optional, for parity tests):**
   ```bash
   cd python
   python create_test_zarr.py
   ```

### Running Tests

1. **All tests (without Zarr feature):**
   ```bash
   cargo test zarr_reader_test
   ```

2. **With Zarr feature enabled:**
   ```bash
   cargo test --features zarr zarr_reader_test
   ```

3. **Including ignored remote tests:**
   ```bash
   cargo test --features zarr zarr_reader_test -- --ignored
   ```

4. **Specific test modules:**
   ```bash
   # Unit tests only
   cargo test --features zarr zarr_reader_test::unit_tests
   
   # Integration tests only
   cargo test --features zarr zarr_reader_test::integration_tests
   
   # Parity tests only
   cargo test --features zarr zarr_reader_test::parity_tests
   ```

## Test Data

### Synthetic Zarr Stores
The tests automatically create minimal synthetic Zarr stores containing:
- **Dimensions**: time (3), level (5), lat (10), lon (12)
- **Coordinates**: time, level, lat, lon with realistic values
- **Variables**: u, v, temp with linear patterns for verification
- **Attributes**: Units, long names, and CF-compliant metadata
- **Compression**: Gzip compression for realistic file structure

### Known Test Patterns
- **U-component wind**: Linear values starting at 10.0, incrementing by 0.1
- **V-component wind**: Linear values starting at 5.0, incrementing by 0.05  
- **Temperature**: Linear values starting at 280.0, incrementing by 0.01
- **Coordinates**: Realistic geographic ranges (30-40°N, 120-110°W)

## Mock Backends

The tests include mock cloud backends that simulate:
- **S3**: `s3://test-bucket/data.zarr`
- **GCS**: `gs://test-bucket/data.zarr`
- **HTTP**: `https://example.com/data.zarr`

These mocks allow testing cloud functionality without requiring actual cloud resources or network access.

## CI Configuration

### Fast CI Tests (enabled by default)
- Unit tests with/without zarr feature
- Integration tests with synthetic data
- Mock cloud backend tests
- Performance tests with small datasets

### Slow CI Tests (ignored by default)
- Remote data access tests
- Large dataset performance tests
- Network-dependent operations

### Test Expectations

1. **Without zarr feature**: All operations should return `UnsupportedOperation` errors
2. **With zarr feature**: Current stub implementation may not fully work but should not panic
3. **Future implementation**: Tests are designed to verify full Zarr functionality once implemented

## Debugging Tests

### Viewing test output:
```bash
cargo test --features zarr zarr_reader_test -- --nocapture
```

### Running individual tests:
```bash
cargo test --features zarr test_zarr_reader_creation -- --exact
```

### Test with debug logging:
```bash
RUST_LOG=debug cargo test --features zarr zarr_reader_test
```

## Adding New Tests

When adding new Zarr functionality:

1. Add unit tests for basic functionality
2. Add integration tests with synthetic data
3. Add parity tests comparing with NetCDF where applicable
4. Add performance tests for large-scale operations
5. Update mock backends if new cloud features are added

## Expected Test Results

### Current Status (with stub implementation)
- ✅ Unit tests: Pass (basic functionality)
- ⚠️ Integration tests: May fail (full reading not implemented)
- ⚠️ Parity tests: May fail (format conversion needed)
- ✅ Mock tests: Pass (testing framework)
- ✅ Performance tests: Pass (metadata access)

### Target Status (with full implementation)
- ✅ All tests should pass
- ✅ Remote data access should work
- ✅ Performance benchmarks should meet targets
