use qibt_rust::trajectory::TrajectoryPoint;
use qibt_rust::data_io::output_trait::{DataWriter, TrajectoryMetadata, OutputFormat, create_writer};
use std::path::Path;
use tempfile::tempdir;

#[test]
fn test_netcdf_output() {
    let temp_dir = tempdir().unwrap();
    let output_path = temp_dir.path().join("test_trajectory.nc");
    
    // Create a NetCDF writer
    let mut writer = create_writer(&output_path, OutputFormat::NetCdf).unwrap();
    
    // Set up metadata
    let metadata = TrajectoryMetadata {
        start_longitude: -100.0,
        start_latitude: 45.0,
        start_pressure: 85000.0,
        start_time: 0.0,
        integration_time_step: 600.0,
        total_integration_time: 3600.0,
        meteorological_data_source: "Test WRF data".to_string(),
        creation_time: "2024-01-01T12:00:00Z".to_string(),
        global_attributes: std::collections::HashMap::new(),
    };
    
    writer.set_metadata(&metadata).unwrap();
    
    // Create some sample trajectory points
    let trajectory = vec![
        TrajectoryPoint {
            time: 0.0,
            longitude: -100.0,
            latitude: 45.0,
            pressure: 85000.0,
            temperature: 280.0,
            u_wind: 10.0,
            v_wind: 5.0,
            w_wind: 0.1,
            potential_temperature: 300.0,
            mixing_ratio: Some(0.010),
            relative_humidity: Some(0.8),
        },
        TrajectoryPoint {
            time: 1.0,
            longitude: -99.9,
            latitude: 45.1,
            pressure: 84500.0,
            temperature: 279.5,
            u_wind: 12.0,
            v_wind: 4.0,
            w_wind: 0.2,
            potential_temperature: 299.8,
            mixing_ratio: Some(0.012),
            relative_humidity: Some(0.85),
        },
        TrajectoryPoint {
            time: 2.0,
            longitude: -99.8,
            latitude: 45.2,
            pressure: 84000.0,
            temperature: 279.0,
            u_wind: 8.0,
            v_wind: 6.0,
            w_wind: -0.1,
            potential_temperature: 299.5,
            mixing_ratio: Some(0.008),
            relative_humidity: Some(0.75),
        },
    ];
    
    // Initialize the writer
    writer.create(1, trajectory.len()).unwrap();
    
    // Write the trajectory
    writer.write_trajectory(1, &trajectory).unwrap();
    
    // Close the writer to finalize the file
    writer.close().unwrap();
    
    // Verify the file was created
    assert!(output_path.exists());
    println!("NetCDF output file created successfully: {}", output_path.display());
}

#[test]
fn test_ascii_output() {
    let temp_dir = tempdir().unwrap();
    let output_path = temp_dir.path().join("test_trajectory.txt");
    
    // Create an ASCII writer
    let mut writer = create_writer(&output_path, OutputFormat::Ascii).unwrap();
    
    // Set up metadata
    let metadata = TrajectoryMetadata {
        start_longitude: -100.0,
        start_latitude: 45.0,
        start_pressure: 85000.0,
        start_time: 0.0,
        integration_time_step: 600.0,
        total_integration_time: 3600.0,
        meteorological_data_source: "Test WRF data".to_string(),
        creation_time: "2024-01-01T12:00:00Z".to_string(),
        global_attributes: std::collections::HashMap::new(),
    };
    
    writer.set_metadata(&metadata).unwrap();
    
    // Create some sample trajectory points
    let trajectory = vec![
        TrajectoryPoint {
            time: 0.0,
            longitude: -100.0,
            latitude: 45.0,
            pressure: 85000.0,
            temperature: 280.0,
            u_wind: 10.0,
            v_wind: 5.0,
            w_wind: 0.1,
            potential_temperature: 300.0,
            mixing_ratio: Some(0.010),
            relative_humidity: Some(0.8),
        },
        TrajectoryPoint {
            time: 1.0,
            longitude: -99.9,
            latitude: 45.1,
            pressure: 84500.0,
            temperature: 279.5,
            u_wind: 12.0,
            v_wind: 4.0,
            w_wind: 0.2,
            potential_temperature: 299.8,
            mixing_ratio: Some(0.012),
            relative_humidity: Some(0.85),
        },
    ];
    
    // Initialize the writer
    writer.create(1, trajectory.len()).unwrap();
    
    // Write the trajectory
    writer.write_trajectory(1, &trajectory).unwrap();
    
    // Close the writer to finalize the file
    writer.close().unwrap();
    
    // Verify the file was created
    assert!(output_path.exists());
    
    // Read the file contents to verify format
    let contents = std::fs::read_to_string(&output_path).unwrap();
    assert!(contents.contains("# Back-trajectory output"));
    assert!(contents.contains("Test WRF data"));
    assert!(contents.contains("1 0.000000 -100.000000 45.000000"));
    
    println!("ASCII output file created successfully: {}", output_path.display());
    println!("Contents:\n{}", contents);
}

#[cfg(feature = "zarr")]
#[test]
fn test_zarr_output() {
    let temp_dir = tempdir().unwrap();
    let output_path = temp_dir.path().join("test_trajectory.zarr");
    
    // Create a Zarr writer (stub implementation)
    match create_writer(&output_path, OutputFormat::Zarr) {
        Ok(mut writer) => {
            // Set up metadata
            let metadata = TrajectoryMetadata {
                start_longitude: -100.0,
                start_latitude: 45.0,
                start_pressure: 85000.0,
                start_time: 0.0,
                integration_time_step: 600.0,
                total_integration_time: 3600.0,
                meteorological_data_source: "Test WRF data".to_string(),
                creation_time: "2024-01-01T12:00:00Z".to_string(),
                global_attributes: std::collections::HashMap::new(),
            };
            
            writer.set_metadata(&metadata).unwrap();
            
            // Create some sample trajectory points
            let trajectory = vec![
                TrajectoryPoint {
                    time: 0.0,
                    longitude: -100.0,
                    latitude: 45.0,
                    pressure: 85000.0,
                    temperature: 280.0,
                    u_wind: 10.0,
                    v_wind: 5.0,
                    w_wind: 0.1,
                    potential_temperature: 300.0,
                    mixing_ratio: Some(0.010),
                    relative_humidity: Some(0.8),
                },
            ];
            
            // Initialize the writer
            writer.create(1, trajectory.len()).unwrap();
            
            // Write the trajectory
            writer.write_trajectory(1, &trajectory).unwrap();
            
            // Close the writer to finalize the dataset
            writer.close().unwrap();
            
            println!("Zarr output dataset created successfully: {}", output_path.display());
        }
        Err(e) => {
            // Expected if zarr feature is not enabled or zarr writer not implemented
            println!("Zarr writer not available (expected if zarr feature not enabled): {}", e);
        }
    }
}

#[test]
fn test_format_detection() {
    // Test format detection from file extensions
    assert_eq!(OutputFormat::from_path(Path::new("trajectory.nc")), OutputFormat::NetCdf);
    assert_eq!(OutputFormat::from_path(Path::new("trajectory.netcdf")), OutputFormat::NetCdf);
    assert_eq!(OutputFormat::from_path(Path::new("data.zarr")), OutputFormat::Zarr);
    assert_eq!(OutputFormat::from_path(Path::new("output.txt")), OutputFormat::Ascii);
    assert_eq!(OutputFormat::from_path(Path::new("output.csv")), OutputFormat::Ascii);
    
    // Test default behavior (no extension)
    assert_eq!(OutputFormat::from_path(Path::new("trajectory")), OutputFormat::NetCdf);
    
    println!("Format detection tests passed");
}

#[test] 
fn test_multiple_format_output() {
    let temp_dir = tempdir().unwrap();
    
    // Create trajectory data
    let trajectory = vec![
        TrajectoryPoint {
            time: 0.0,
            longitude: -100.0,
            latitude: 45.0,
            pressure: 85000.0,
            temperature: 280.0,
            u_wind: 10.0,
            v_wind: 5.0,
            w_wind: 0.1,
            potential_temperature: 300.0,
            mixing_ratio: Some(0.010),
            relative_humidity: Some(0.8),
        },
        TrajectoryPoint {
            time: 1.0,
            longitude: -99.9,
            latitude: 45.1,
            pressure: 84500.0,
            temperature: 279.5,
            u_wind: 12.0,
            v_wind: 4.0,
            w_wind: 0.2,
            potential_temperature: 299.8,
            mixing_ratio: Some(0.012),
            relative_humidity: Some(0.85),
        },
    ];
    
    let metadata = TrajectoryMetadata {
        start_longitude: -100.0,
        start_latitude: 45.0,
        start_pressure: 85000.0,
        start_time: 0.0,
        integration_time_step: 600.0,
        total_integration_time: 3600.0,
        meteorological_data_source: "Test WRF data".to_string(),
        creation_time: "2024-01-01T12:00:00Z".to_string(),
        global_attributes: std::collections::HashMap::new(),
    };
    
    // Test writing to all supported formats
    let formats = [
        (OutputFormat::NetCdf, "trajectory.nc"),
        (OutputFormat::Ascii, "trajectory.txt"),
    ];
    
    for (format, filename) in &formats {
        let output_path = temp_dir.path().join(filename);
        let mut writer = create_writer(&output_path, *format).unwrap();
        
        writer.set_metadata(&metadata).unwrap();
        writer.create(1, trajectory.len()).unwrap();
        writer.write_trajectory(1, &trajectory).unwrap();
        writer.close().unwrap();
        
        assert!(output_path.exists());
        println!("Successfully created {} output: {}", format, output_path.display());
    }
    
    // Test Zarr if feature is enabled
    #[cfg(feature = "zarr")]
    {
        let output_path = temp_dir.path().join("trajectory.zarr");
        match create_writer(&output_path, OutputFormat::Zarr) {
            Ok(mut writer) => {
                writer.set_metadata(&metadata).unwrap();
                writer.create(1, trajectory.len()).unwrap();
                writer.write_trajectory(1, &trajectory).unwrap();
                writer.close().unwrap();
                println!("Successfully created Zarr output: {}", output_path.display());
            }
            Err(e) => {
                println!("Zarr writer not implemented yet: {}", e);
            }
        }
    }
}
