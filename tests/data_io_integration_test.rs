use ndarray::{Array1, Array2, Array4};
use qibt_rust::data_io::{MeteoData, MeteoField, MeteoFieldArray, MeteoGrid};

#[test]
fn test_meteo_field_creation() {
    let field = MeteoField {
        name: "temperature".to_string(),
        data: vec![vec![vec![vec![273.15, 274.15], vec![275.15, 276.15]]]],
        missing_value: -9999.0,
        units: "K".to_string(),
    };

    assert_eq!(field.name, "temperature");
    assert_eq!(field.units, "K");
    assert_eq!(field.missing_value, -9999.0);
}

#[test]
fn test_meteo_field_array_creation() {
    let data = Array4::zeros((10, 15, 20, 24)); // [j, i, k, t]
    let field = MeteoFieldArray {
        name: "U".to_string(),
        data,
        missing_value: -9999.0,
        units: "m/s".to_string(),
    };

    assert_eq!(field.name, "U");
    assert_eq!(field.units, "m/s");
    assert_eq!(field.missing_value, -9999.0);
    assert_eq!(field.data.shape(), &[10, 15, 20, 24]);
}

#[test]
fn test_meteo_grid_creation() {
    let grid = MeteoGrid {
        longitudes: vec![-180.0, -179.0, -178.0],
        latitudes: vec![90.0, 89.0, 88.0],
        levels: vec![100000.0, 95000.0, 90000.0],
        times: vec![0.0, 1.0, 2.0],
        time_reference: "hours since 2023-01-01 00:00:00".to_string(),
    };

    assert_eq!(grid.longitudes.len(), 3);
    assert_eq!(grid.latitudes.len(), 3);
    assert_eq!(grid.levels.len(), 3);
    assert_eq!(grid.times.len(), 3);
}

#[test]
fn test_meteo_data_operations() {
    let mut data = MeteoData::new();

    // Test that new MeteoData is empty
    assert_eq!(data.fields.len(), 0);
    assert_eq!(data.grid.longitudes.len(), 0);

    // Add a field
    let field = MeteoField {
        name: "u_wind".to_string(),
        data: vec![],
        missing_value: -9999.0,
        units: "m/s".to_string(),
    };

    data.add_field(field);
    assert_eq!(data.fields.len(), 1);

    // Test field retrieval
    let retrieved = data.get_field("u_wind");
    assert!(retrieved.is_some());
    assert_eq!(retrieved.unwrap().name, "u_wind");

    let missing = data.get_field("nonexistent");
    assert!(missing.is_none());
}

#[test]
fn test_array_dimension_layout() {
    // Test that our array layout matches expectations for interpolation [j, i, k, t]
    let nj = 50; // south-north (j direction)
    let ni = 60; // west-east (i direction) 
    let nk = 30; // vertical levels (k direction)
    let nt = 24; // time steps (t direction)

    let array = Array4::<f32>::zeros((nj, ni, nk, nt));

    assert_eq!(array.shape()[0], nj); // j dimension
    assert_eq!(array.shape()[1], ni); // i dimension  
    assert_eq!(array.shape()[2], nk); // k dimension
    assert_eq!(array.shape()[3], nt); // t dimension

    // Test that we can access elements in the expected order
    let j_idx = 25;
    let i_idx = 30;
    let k_idx = 15;
    let t_idx = 12;

    // This should not panic and should access the correct element
    let _value = array[[j_idx, i_idx, k_idx, t_idx]];
}

#[test]
fn test_coordinate_array_shapes() {
    // Test typical WRF coordinate array shapes
    let ni = 100;
    let nj = 80;
    let nk = 35;

    // Lat/Lon are typically 2D (ni, nj)
    let xlat = Array2::<f32>::zeros((nj, ni));
    let xlong = Array2::<f32>::zeros((nj, ni));

    assert_eq!(xlat.shape(), &[nj, ni]);
    assert_eq!(xlong.shape(), &[nj, ni]);

    // Pressure levels are typically 1D
    let levels = Array1::<f32>::zeros(nk);
    assert_eq!(levels.shape(), &[nk]);
}

#[test]
fn test_boundary_trimming_simulation() {
    // Simulate what boundary trimming should do
    let original_shape = (52, 62, 30, 24); // Original dimensions
    let boundary_trim = 1;

    let _original_array = Array4::<f32>::ones(original_shape);

    // Simulate boundary trimming: remove 1 element from each side of j and i dimensions
    let expected_shape = (
        original_shape.0 - 2 * boundary_trim, // j dimension: 52 - 2 = 50
        original_shape.1 - 2 * boundary_trim, // i dimension: 62 - 2 = 60
        original_shape.2,                     // k dimension: unchanged
        original_shape.3,                     // t dimension: unchanged
    );

    // This simulates what the reader should produce after trimming
    let trimmed_shape = (
        expected_shape.0,
        expected_shape.1,
        expected_shape.2,
        expected_shape.3,
    );
    let trimmed_array = Array4::<f32>::ones(trimmed_shape);

    assert_eq!(trimmed_array.shape(), &[50, 60, 30, 24]);
}

#[test]
fn test_standard_variable_names() {
    // Test that we're using the correct WRF variable names
    let standard_vars = ["U", "V", "W", "T", "QVAPOR", "RAIN"];

    for var in &standard_vars {
        assert!(!var.is_empty());
        assert!(
            var.chars()
                .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit())
        );
    }

    assert_eq!(standard_vars.len(), 6);
}
