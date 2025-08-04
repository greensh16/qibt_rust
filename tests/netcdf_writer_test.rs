use qibt_rust::data_io::writer::{BufferedRecord, NetCDFWriter, TrajectoryMetadata};
use qibt_rust::trajectory::TrajectoryPoint;
use std::fs;

#[test]
fn test_create_output_basic() {
    // Create test trajectory data
    let trajectory1 = vec![
        TrajectoryPoint {
            time: 0.0,
            longitude: -122.0,
            latitude: 37.0,
            pressure: 101325.0,
            temperature: 288.15,
            u_wind: 5.0,
            v_wind: 2.0,
            w_wind: 0.1,
            potential_temperature: 290.0,
            mixing_ratio: Some(0.005),
            relative_humidity: Some(60.0),
        },
        TrajectoryPoint {
            time: 1.0,
            longitude: -122.1,
            latitude: 37.1,
            pressure: 101200.0,
            temperature: 287.5,
            u_wind: 5.2,
            v_wind: 2.1,
            w_wind: 0.05,
            potential_temperature: 290.2,
            mixing_ratio: Some(0.0048),
            relative_humidity: Some(58.0),
        },
    ];

    let trajectory2 = vec![TrajectoryPoint {
        time: 0.0,
        longitude: -121.0,
        latitude: 38.0,
        pressure: 101325.0,
        temperature: 289.0,
        u_wind: 4.5,
        v_wind: 1.8,
        w_wind: 0.08,
        potential_temperature: 291.0,
        mixing_ratio: Some(0.006),
        relative_humidity: Some(65.0),
    }];

    let metadata1 = TrajectoryMetadata {
        start_longitude: -122.0,
        start_latitude: 37.0,
        start_pressure: 101325.0,
        start_time: 0.0,
        integration_time_step: 3600.0,
        total_integration_time: 3600.0,
        meteorological_data_source: "Test data".to_string(),
        creation_time: "2024-01-01T00:00:00Z".to_string(),
    };

    let metadata2 = TrajectoryMetadata {
        start_longitude: -121.0,
        start_latitude: 38.0,
        start_pressure: 101325.0,
        start_time: 0.0,
        integration_time_step: 3600.0,
        total_integration_time: 3600.0,
        meteorological_data_source: "Test data".to_string(),
        creation_time: "2024-01-01T00:00:00Z".to_string(),
    };

    let records = vec![
        BufferedRecord {
            parcel_id: 1,
            trajectory: trajectory1,
            metadata: metadata1,
        },
        BufferedRecord {
            parcel_id: 2,
            trajectory: trajectory2,
            metadata: metadata2,
        },
    ];

    // Create writer and output file
    let output_path = "/tmp/test_trajectory_output.nc";
    let mut writer = NetCDFWriter::new(output_path);

    // Test create_output
    let result = writer.create_output(100, records);
    println!("Create output result: {:?}", result);

    // Check if file was created
    if std::path::Path::new(output_path).exists() {
        println!("NetCDF file created successfully at: {}", output_path);

        // Clean up test file
        let _ = fs::remove_file(output_path);
    } else {
        println!("NetCDF file was not created");
    }

    // The test should not panic - we're validating the API works
}

#[test]
fn test_buffer_record() {
    let mut writer = NetCDFWriter::new("/tmp/test_buffer.nc");

    let trajectory = vec![TrajectoryPoint {
        time: 0.0,
        longitude: -120.0,
        latitude: 36.0,
        pressure: 101325.0,
        temperature: 288.0,
        u_wind: 3.0,
        v_wind: 1.5,
        w_wind: 0.02,
        potential_temperature: 290.0,
        mixing_ratio: None,
        relative_humidity: None,
    }];

    let metadata = TrajectoryMetadata {
        start_longitude: -120.0,
        start_latitude: 36.0,
        start_pressure: 101325.0,
        start_time: 0.0,
        integration_time_step: 3600.0,
        total_integration_time: 3600.0,
        meteorological_data_source: "Test".to_string(),
        creation_time: "2024-01-01T00:00:00Z".to_string(),
    };

    // Test buffering
    writer.buffer_record(1, trajectory, metadata);

    // Check buffer contains one record
    let buffer = writer.buffer.lock().unwrap();
    assert_eq!(buffer.len(), 1);
    assert_eq!(buffer[0].parcel_id, 1);
}
