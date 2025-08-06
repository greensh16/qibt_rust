use qibt_rust::io::{Dataset, DataReader};
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <zarr_path>", args[0]);
        std::process::exit(1);
    }
    
    let zarr_path = &args[1];  
    println!("Testing Zarr store at: {}", zarr_path);
    
    // Try to open the dataset
    println!("Attempting to open Zarr dataset...");
    match Dataset::open(zarr_path) {
        Ok(dataset) => {
            println!("✓ Successfully opened Zarr dataset");
            
            // List variables
            match dataset.list_variables() {
                Ok(vars) => {
                    println!("Variables found: {:?}", vars);
                },
                Err(e) => {
                    println!("✗ Error listing variables: {}", e);
                }
            }
            
            // List dimensions
            match dataset.list_dimensions() {
                Ok(dims) => {
                    println!("Dimensions found: {:?}", dims);
                },
                Err(e) => {
                    println!("✗ Error listing dimensions: {}", e);
                }
            }
            
            // Get metadata
            match dataset.get_metadata() {
                Ok(metadata) => {
                    println!("Metadata successfully retrieved:");
                    println!("  Dimensions: {}", metadata.dimensions.len());
                    println!("  Variables: {}", metadata.variables.len());
                    println!("  Global attributes: {}", metadata.global_attributes.len());
                },
                Err(e) => {
                    println!("✗ Error getting metadata: {}", e);
                }
            }
        },
        Err(e) => {
            println!("✗ Failed to open Zarr dataset: {}", e);
        }
    }
    
    Ok(())
}
