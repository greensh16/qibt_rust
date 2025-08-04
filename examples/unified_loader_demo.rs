/// Demonstration of the unified loader functionality
/// 
/// This example shows how the Dataset::open() method can automatically detect
/// and handle different file formats (NetCDF and Zarr) using the unified interface.

use qibt_rust::io::{Dataset, DataReader, DataReaderError, is_netcdf_format, is_zarr_format};
use std::path::Path;

fn main() -> Result<(), DataReaderError> {
    println!("=== Unified Data Loader Demonstration ===\n");
    
    // Example paths - these would be real files in practice
    let example_paths = vec![
        "/path/to/weather_data.nc",     // NetCDF file
        "/path/to/climate_data.zarr",   // Zarr dataset with .zarr extension
        "/path/to/zarr_directory",      // Zarr dataset as directory
        "/path/to/model_output.nc4",    // NetCDF-4 file
    ];
    
    println!("1. Format Detection Examples:");
    println!("============================\n");
    
    // Demonstrate format detection logic
    for path_str in &example_paths {
        let path = Path::new(path_str);
        println!("Checking: {}", path_str);
        
        // Note: These will return errors for non-existent files, but demonstrate the logic
        match is_netcdf_format(path) {
            Ok(true) => println!("  → Detected as NetCDF format"),
            Ok(false) => {
                match is_zarr_format(path) {
                    Ok(true) => println!("  → Detected as Zarr format"),
                    Ok(false) => println!("  → Format not detected (would fall back to extensions)"),
                    Err(_) => println!("  → File not found (demonstration only)"),
                }
            },
            Err(_) => println!("  → File not found (demonstration only)"),
        }
        println!();
    }
    
    println!("2. Unified Dataset Opening:");
    println!("===========================\n");
    
    // Demonstrate the unified Dataset::open() interface
    for path_str in &example_paths {
        println!("Attempting to open: {}", path_str);
        
        // This is how users would open datasets with automatic format detection
        match Dataset::open(path_str) {
            Ok(dataset) => {
                println!("  ✓ Successfully opened dataset");
                println!("  → Reader is open: {}", dataset.is_open());
                
                // Example of using the dataset
                match dataset.list_variables() {
                    Ok(vars) => println!("  → Variables: {:?}", vars),
                    Err(e) => println!("  → Could not list variables: {}", e),
                }
            },
            Err(e) => {
                println!("  ✗ Failed to open: {}", e);
                // In real usage, this would be due to actual file access issues
                // or unsupported formats, not just missing files
            }
        }
        println!();
    }
    
    println!("3. Implementation Details:");
    println!("=========================\n");
    
    println!("The unified loader implements the following detection strategy:");
    println!("  1. First tries NetCDF magic bytes detection:");
    println!("     - Classic NetCDF: 'CDF\\001' or 'CDF\\002'");
    println!("     - NetCDF-4 (HDF5): '\\211HDF\\r\\n\\032\\n'");
    println!("  2. If NetCDF detection fails, checks for Zarr format:");
    println!("     - '.zarr' suffix in path name");
    println!("     - Directory containing '.zgroup' file (Zarr group)");
    println!("     - Directory containing '.zarray' file (Zarr array)"); 
    println!("     - JSON file with Zarr metadata fields");
    println!("  3. Falls back to extension-based detection for compatibility");
    println!();
    
    println!("Key Features:");
    println!("  ✓ Preserves existing NetCDF behavior unchanged");
    println!("  ✓ Returns appropriate Box<dyn DataReader> for each format");
    println!("  ✓ Provides unified Dataset interface for both formats");
    println!("  ✓ Supports both file paths and URLs (when readers support them)");
    println!("  ✓ Graceful error handling with descriptive messages");
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_demonstration_logic() {
        // This test verifies that the demonstration code logic is sound
        // Even though files don't exist, the error handling should work
        
        let result = Dataset::open("/nonexistent/file.nc");
        assert!(result.is_err());
        
        // The error should be descriptive
        let error_msg = result.unwrap_err().to_string();
        assert!(error_msg.contains("Unable to detect format") || error_msg.contains("File not found"));
    }
}
