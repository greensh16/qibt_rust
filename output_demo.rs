use qibt_rust::data_io::output_trait::{TrajectoryMetadata, OutputFormat, create_writer};
use qibt_rust::trajectory::Trajectory;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Creating sample trajectory output files...");
    
    // Create sample trajectory data
    let mut trajectory = Trajectory::new();
    trajectory.time_hours = vec![0.0, 1.0, 2.0, 3.0];
    trajectory.longitude = vec![-100.0, -99.9, -99.8, -99.7];
    trajectory.latitude = vec![45.0, 45.1, 45.2, 45.3];
    trajectory.pressure = vec![85000.0, 84500.0, 84000.0, 83500.0];
    trajectory.temperature = vec![280.0, 279.5, 279.0, 278.5];
    trajectory.u_wind = vec![10.0, 12.0, 11.0, 9.0];
    trajectory.v_wind = vec![5.0, 4.0, 6.0, 7.0];
    trajectory.w_wind = vec![0.1, 0.2, -0.1, 0.0];
    
    let trajectories = vec![trajectory];
    
    // Create metadata
    let metadata = TrajectoryMetadata {
        creation_time: "2024-08-05T11:33:43Z".to_string(),
        data_source: "Demo WRF data".to_string(),
        start_location: (-100.0, 45.0),
        start_pressure: 85000.0,
    };
    
    // Create output files in different formats
    let output_dir = "output_samples";
    std::fs::create_dir_all(output_dir)?;
    
    // 1. NetCDF output
    let netcdf_path = format!("{}/sample_trajectory.nc", output_dir);
    println!("Creating NetCDF file: {}", netcdf_path);
    let mut netcdf_writer = create_writer(OutputFormat::NetCDF, &netcdf_path)?;
    netcdf_writer.write_trajectories(&trajectories, &metadata)?;
    netcdf_writer.close()?;
    println!("✓ NetCDF file created successfully");
    
    // 2. ASCII output
    let ascii_path = format!("{}/sample_trajectory.txt", output_dir);
    println!("Creating ASCII file: {}", ascii_path);
    let mut ascii_writer = create_writer(OutputFormat::ASCII, &ascii_path)?;
    ascii_writer.write_trajectories(&trajectories, &metadata)?;
    ascii_writer.close()?;
    println!("✓ ASCII file created successfully");
    
    // 3. Try Zarr output (if feature is enabled)
    #[cfg(feature = "zarr")] 
    {
        let zarr_path = format!("{}/sample_trajectory.zarr", output_dir);
        println!("Creating Zarr store: {}", zarr_path);
        let mut zarr_writer = create_writer(OutputFormat::Zarr, &zarr_path)?;
        zarr_writer.write_trajectories(&trajectories, &metadata)?;
        zarr_writer.close()?;
        println!("✓ Zarr store created successfully");
    }
    
    #[cfg(not(feature = "zarr"))]
    {
        println!("⚠ Zarr output not available (zarr feature not enabled)");
    }
    
    println!("\nOutput files created in '{}':", output_dir);
    println!("  - sample_trajectory.nc (NetCDF format)");
    println!("  - sample_trajectory.txt (ASCII format)");
    #[cfg(feature = "zarr")]
    println!("  - sample_trajectory.zarr/ (Zarr format)");
    
    // Display content of ASCII file
    println!("\n--- ASCII File Content ---");
    let ascii_content = std::fs::read_to_string(&ascii_path)?;
    println!("{}", ascii_content);
    
    Ok(())
}
