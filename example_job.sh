#!/bin/bash
#
# Example job script for QIBT Rust
# 
# This script demonstrates how to run the compiled QIBT Rust binary
# with typical parameters for trajectory computation.
#
# Make sure to:
# 1. Build the project first: cargo build --release
# 2. Adjust input/output paths as needed
# 3. Create output directory if it doesn't exist
#

# Configuration
INPUT_FILE="test_data/wrfout_d01_2023-07-31_12:00:00.nc"
OUTPUT_FILE="results/output_trajectory.nc"
NUM_PARCELS=100
NUM_THREADS=8
START_LAT=45.0
START_LON=-100.0
TRAJECTORY_LENGTH=120  # hours

# Create output directory if it doesn't exist
mkdir -p $(dirname "$OUTPUT_FILE")

# Check if input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' does not exist."
    echo "Please update the INPUT_FILE variable or provide the correct path."
    exit 1
fi

# Check if binary exists
BINARY="./target/release/qibt_rust"
if [[ ! -f "$BINARY" ]]; then
    echo "Error: Binary '$BINARY' does not exist."
    echo "Please build the project first: cargo build --release"
    exit 1
fi

echo "Starting QIBT Rust trajectory computation..."
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo "Parcels: $NUM_PARCELS"
echo "Threads: $NUM_THREADS"
echo "Starting position: ($START_LAT, $START_LON)"
echo "Trajectory length: $TRAJECTORY_LENGTH hours"
echo ""

# Run the trajectory computation
$BINARY run \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_FILE" \
    --parcels $NUM_PARCELS \
    --threads $NUM_THREADS \
    --start-lat $START_LAT \
    --start-lon $START_LON \
    --length $TRAJECTORY_LENGTH

# Check if the command was successful
if [[ $? -eq 0 ]]; then
    echo ""
    echo "Trajectory computation completed successfully!"
    echo "Output saved to: $OUTPUT_FILE"
else
    echo ""
    echo "Error: Trajectory computation failed."
    exit 1
fi
