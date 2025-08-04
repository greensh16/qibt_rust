# QIBT Rust

Quasi-Isentropic Back-Trajectory Analysis Tool in Rust.

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

### Example:

```bash
./target/release/qibt_rust run --input <input_file> --output <output_file>
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