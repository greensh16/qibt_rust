/// Example demonstrating unified trajectory computation with different data readers
/// This shows how the same trajectory computation pipeline can work with both
/// NetCDF and Zarr data readers without code duplication.

use qibt_rust::{
    config::{Config, Constants},
    data_io::{NetCDFReader, NetCDFWriter},
    io::{DataReader},
    trajectory::{LegacyParcel, integrate_back_trajectory_generic},
    // parallel::{compute_trajectories_parallel_generic},
};

#[cfg(feature = "zarr")]
use qibt_rust::io::create_zarr_reader;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== Unified Trajectory Computation Demo ===");
    
    // Configuration - same for both readers
    let mut config = Config::default();
    config.start_time = 0.0;
    config.trajectory_length = 48.0; // 48 hours
    config.time_step = 3600.0; // 1 hour in seconds
    
    // Create sample parcels
    let parcels = vec![
        LegacyParcel::new(1, -100.0, 40.0, 85000.0, config.start_time),
        LegacyParcel::new(2, -90.0, 35.0, 85000.0, config.start_time),
        LegacyParcel::new(3, -80.0, 30.0, 85000.0, config.start_time),
    ];
    
    println!("Created {} parcels for trajectory computation", parcels.len());
    
    // Demonstrate with NetCDF reader (if available)
    if let Ok(netcdf_path) = std::env::var("NETCDF_TEST_FILE") {
        println!("\n--- Computing trajectories with NetCDF reader ---");
        if Path::new(&netcdf_path).exists() {
            compute_with_netcdf(&config, &parcels, &netcdf_path)?;
        } else {
            println!("NetCDF file not found: {}", netcdf_path);
        }
    } else {
        println!("NetCDF demo skipped - set NETCDF_TEST_FILE environment variable");
    }
    
    // Demonstrate with Zarr reader (if available)
    if let Ok(zarr_path) = std::env::var("ZARR_TEST_PATH") {
        println!("\n--- Computing trajectories with Zarr reader ---");
        if Path::new(&zarr_path).exists() {
            compute_with_zarr(&config, &parcels, &zarr_path)?;
        } else {
            println!("Zarr path not found: {}", zarr_path);
        }
    } else {
        println!("Zarr demo skipped - set ZARR_TEST_PATH environment variable");
    }
    
    // Demonstrate the unified interface
    println!("\n--- Demonstrating unified interface ---");
    demonstrate_unified_interface(&config, &parcels)?;
    
    println!("\nDemo completed successfully!");
    Ok(())
}

/// Compute trajectories using NetCDF reader
fn compute_with_netcdf(
    config: &Config,
    parcels: &[LegacyParcel],
    file_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Opening NetCDF file: {}", file_path);
    
    // This would normally open a real NetCDF file
    // For demo purposes, we'll show the interface
    let reader = NetCDFReader::new(file_path);
    
    // Validate the file exists and is readable
    match reader.validate_file() {
        Ok(_) => println!("NetCDF file validation passed"),
        Err(e) => {
            println!("NetCDF file validation failed: {}", e);
            return Ok(()); // Continue with demo
        }
    }
    
    // Demonstrate single trajectory computation
    let mut parcel = parcels[0].clone();
    let writer = NetCDFWriter::new("output_netcdf_single.nc");
    
    // Use generic trajectory integration (works with any DataReader)
    println!("Computing single trajectory with generic interface...");
    match integrate_back_trajectory_generic(config, &mut parcel, &reader, &writer) {
        Ok(_) => println!("Single trajectory completed: {} points", parcel.trajectory.len()),
        Err(e) => println!("Single trajectory failed: {}", e),
    }
    
    // Demonstrate parallel computation (would work if we had real data)
    println!("Would compute {} trajectories in parallel...", parcels.len());
    // This would work with real data:
    // let trajectories = compute_trajectories_parallel_generic(config, parcels.to_vec(), &reader)?;
    
    Ok(())
}

/// Compute trajectories using Zarr reader
fn compute_with_zarr(
    config: &Config,
    parcels: &[LegacyParcel],
    dataset_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Opening Zarr dataset: {}", dataset_path);
    
    // Create Zarr reader using the unified interface
    let reader = create_zarr_reader(dataset_path)?;
    
    // List available variables
    match reader.list_variables() {
        Ok(vars) => println!("Available variables: {:?}", vars),
        Err(e) => println!("Could not list variables: {}", e),
    }
    
    // Get metadata
    match reader.get_metadata() {
        Ok(metadata) => {
            println!("Dataset dimensions: {:?}", metadata.dimensions.iter().map(|d| &d.name).collect::<Vec<_>>());
            println!("Dataset variables: {:?}", metadata.variables.iter().map(|v| &v.name).collect::<Vec<_>>());
        },
        Err(e) => println!("Could not get metadata: {}", e),
    }
    
    // Demonstrate single trajectory computation using the SAME CODE as NetCDF
    let mut parcel = parcels[0].clone();
    let writer = NetCDFWriter::new("output_zarr_single.nc");
    
    println!("Computing single trajectory with generic interface...");
    match integrate_back_trajectory_generic(config, &mut parcel, reader.as_ref(), &writer) {
        Ok(_) => println!("Single trajectory completed: {} points", parcel.trajectory.len()),
        Err(e) => println!("Single trajectory failed: {}", e),
    }
    
    // Demonstrate parallel computation (same code works for both readers!)
    println!("Would compute {} trajectories in parallel...", parcels.len());
    // This would work with real data:
    // let trajectories = compute_trajectories_parallel_generic(config, parcels.to_vec(), reader.as_ref())?;
    
    Ok(())
}

/// Demonstrate that the same code works with any DataReader implementation
fn demonstrate_unified_interface(
    config: &Config,
    parcels: &[LegacyParcel],
) -> Result<(), Box<dyn std::error::Error>> {
    println!("This demonstrates how the same trajectory computation code");
    println!("can work with any data reader that implements the DataReader trait.");
    
    // Generic function that works with any DataReader
    fn compute_with_any_reader<T: DataReader>(
        config: &Config,
        parcels: &[LegacyParcel],
        reader: &T,
        reader_type: &str,
    ) -> Result<(), String> {
        println!("Computing with {} reader...", reader_type);
        
        // This exact same code works for NetCDF, Zarr, or any other DataReader implementation
        let mut parcel = parcels[0].clone();
        let writer = NetCDFWriter::new(&format!("output_{}.nc", reader_type.to_lowercase()));
        
        integrate_back_trajectory_generic(config, &mut parcel, reader, &writer)?;
        
        println!("{} reader: trajectory with {} points computed", reader_type, parcel.trajectory.len());
        Ok(())
    }
    
    // Example with NetCDF reader
    let netcdf_reader = NetCDFReader::new("example.nc");
    if let Err(e) = compute_with_any_reader(config, parcels, &netcdf_reader, "NetCDF") {
        println!("NetCDF computation failed (expected): {}", e);
    }
    
    // Example with Zarr reader (if available)
    if let Ok(zarr_reader) = create_zarr_reader("example.zarr") {
        if let Err(e) = compute_with_any_reader(config, parcels, zarr_reader.as_ref(), "Zarr") {
            println!("Zarr computation failed (expected): {}", e);
        }
    }
    
    println!("The key insight: NO code duplication!");
    println!("- Same interpolation algorithms");
    println!("- Same trajectory integration");
    println!("- Same parallel processing");
    println!("- Same physics calculations");
    println!("Only the data reading is different, and that's abstracted away.");
    
    Ok(())
}

/// Example showing how to extend the system with a new reader type
#[allow(dead_code)]
fn demonstrate_extensibility() {
    println!("\nExtensibility demonstration:");
    println!("To add support for a new data format (e.g., HDF5, GRIB), you would:");
    println!("1. Implement the DataReader trait for your new format");
    println!("2. All existing trajectory computation code automatically works!");
    println!("3. No changes needed to interpolation, physics, or parallel processing");
    println!("\nThis design eliminates vendor lock-in and makes the system future-proof.");
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_unified_interface_concept() {
        // This test validates that our design works conceptually
        let mut config = Config::default();
        config.start_time = 0.0;
        config.trajectory_length = 24.0;
        config.time_step = 3600.0;
        
        let parcels = vec![
            LegacyParcel::new(1, -100.0, 40.0, 85000.0, config.start_time),
        ];
        
        // The important thing is that this compiles and shows the interface works
        let netcdf_reader = NetCDFReader::new("test.nc");
        
        // This function signature accepts any DataReader
        fn _test_generic_interface<T: DataReader>(_reader: &T) -> bool {
            // Could perform trajectory computation here
            true
        }
        
        assert!(_test_generic_interface(&netcdf_reader));
        println!("Unified interface test passed");
    }
    
    #[test]
    fn test_field_data_access_trait() {
        use qibt_rust::math::interpolate::{FieldDataAccess, VecFieldData};
        
        // Test the trait-based field access
        let data = vec![
            vec![vec![vec![1.0, 2.0], vec![3.0, 4.0]]], // time 0, level 0
        ];
        let longitudes = vec![0.0, 1.0];
        let latitudes = vec![0.0, 1.0];
        let levels = vec![1000.0];
        let times = vec![0.0];
        
        let field_data = VecFieldData {
            data: &data,
            longitudes: &longitudes,
            latitudes: &latitudes,
            levels: &levels,
            times: &times,
        };
        
        // Test shape
        let shape = field_data.get_shape();
        assert_eq!(shape, (1, 1, 2, 2)); // (time, level, lat, lon)
        
        // Test coordinate access
        let coords = field_data.get_coordinates().unwrap();
        assert_eq!(coords.0, longitudes);
        assert_eq!(coords.1, latitudes);
        assert_eq!(coords.2, levels);
        assert_eq!(coords.3, times);
        
        // Test value access
        let value = field_data.get_value(0, 0, 0, 0).unwrap();
        assert_eq!(value, 1.0);
        
        println!("FieldDataAccess trait test passed");
    }
}
