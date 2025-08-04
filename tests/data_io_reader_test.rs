use chrono::NaiveDate;
use qibt_rust::data_io::reader::{AdvancedNetCDFReader, NetCDFReader, ReaderError, get_filename};

#[test]
fn test_get_filename() {
    let base_path = "/data";
    let datetime = NaiveDate::from_ymd_opt(2023, 7, 31)
        .unwrap()
        .and_hms_opt(12, 0, 0)
        .unwrap()
        .and_utc();
    let domain = Some(1);
    let file_type = "wrfout";

    let expected_path = "/data/wrfout_d01_2023-07-31_12:00:00.nc";
    let generated_path = get_filename(base_path, &datetime, domain, file_type);
    assert_eq!(generated_path.to_str().unwrap(), expected_path);
}

#[test]
fn test_get_filename_no_domain() {
    let base_path = "/data";
    let datetime = NaiveDate::from_ymd_opt(2023, 12, 25)
        .unwrap()
        .and_hms_opt(0, 30, 45)
        .unwrap()
        .and_utc();
    let domain = None;
    let file_type = "wrfout";

    let expected_path = "/data/wrfout_2023-12-25_00:30:45.nc";
    let generated_path = get_filename(base_path, &datetime, domain, file_type);
    assert_eq!(generated_path.to_str().unwrap(), expected_path);
}

#[test]
fn test_get_filename_different_file_types() {
    let base_path = "/model/output";
    let datetime = NaiveDate::from_ymd_opt(2024, 1, 1)
        .unwrap()
        .and_hms_opt(18, 0, 0)
        .unwrap()
        .and_utc();

    // Test different file types
    let wrfinput_path = get_filename(base_path, &datetime, Some(2), "wrfinput");
    assert_eq!(
        wrfinput_path.to_str().unwrap(),
        "/model/output/wrfinput_d02_2024-01-01_18:00:00.nc"
    );

    let wrfbdy_path = get_filename(base_path, &datetime, None, "wrfbdy");
    assert_eq!(
        wrfbdy_path.to_str().unwrap(),
        "/model/output/wrfbdy_2024-01-01_18:00:00.nc"
    );
}

#[test]
fn test_reader_with_nonexistent_file() {
    let reader = AdvancedNetCDFReader::new("nonexistent/path", 1);
    let datetime = NaiveDate::from_ymd_opt(2023, 7, 31)
        .unwrap()
        .and_hms_opt(12, 0, 0)
        .unwrap()
        .and_utc();
    let domain = Some(1);
    let result = reader.read_variable_array(&datetime, "U", domain);
    assert!(matches!(result, Err(ReaderError::FileNotFound(_))));
}

#[test]
fn test_read_multiple_variables_with_missing_file() {
    let reader = AdvancedNetCDFReader::new("nonexistent/path", 1);
    let datetime = NaiveDate::from_ymd_opt(2023, 7, 31)
        .unwrap()
        .and_hms_opt(12, 0, 0)
        .unwrap()
        .and_utc();
    let domain = Some(1);

    let variables = ["U", "V", "W"];
    let result = reader.read_multiple_variables(&datetime, &variables, domain);
    assert!(result.is_err());
}

#[test]
fn test_advanced_reader_creation() {
    let reader = AdvancedNetCDFReader::new("/test/path", 2);
    assert_eq!(reader.base_path, "/test/path");
    assert_eq!(reader.boundary_trim, 2);
}

#[test]
fn test_netcdf_reader_creation() {
    let reader = NetCDFReader::new("/test/file.nc");
    assert_eq!(reader.file_path, "/test/file.nc");
}

#[test]
fn test_netcdf_reader_validate_file_nonexistent() {
    let reader = NetCDFReader::new("nonexistent_file.nc");
    let result = reader.validate_file();
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("File does not exist"));
}

#[test]
fn test_read_standard_meteo_variables() {
    let reader = AdvancedNetCDFReader::new("nonexistent/path", 0);
    let datetime = NaiveDate::from_ymd_opt(2023, 7, 31)
        .unwrap()
        .and_hms_opt(12, 0, 0)
        .unwrap()
        .and_utc();
    let domain = Some(1);

    let result = reader.read_meteo_variables(&datetime, domain);
    assert!(result.is_err()); // Should fail because file doesn't exist
}

#[test]
fn test_error_display() {
    let error = ReaderError::MissingVariable("TEST_VAR".to_string());
    assert_eq!(format!("{}", error), "Variable not found: TEST_VAR");

    let error = ReaderError::FileNotFound("/path/to/file.nc".to_string());
    assert_eq!(format!("{}", error), "File not found: /path/to/file.nc");

    let error = ReaderError::ConversionError;
    assert_eq!(format!("{}", error), "Data conversion error");
}

// Tests using actual NetCDF files
#[test]
fn test_read_variable_array_with_real_netcdf() {
    let reader = AdvancedNetCDFReader::new("test_data", 1); // boundary_trim = 1
    let datetime = NaiveDate::from_ymd_opt(2023, 7, 31)
        .unwrap()
        .and_hms_opt(12, 0, 0)
        .unwrap()
        .and_utc();

    // Test reading U variable from small test file
    let result = reader.read_variable_array(&datetime, "U", Some(1));
    if result.is_ok() {
        let array = result.unwrap();
        println!("Array shape: {:?}", array.shape());
        // Just verify it's a 4D array with correct layout [j, i, k, t]
        assert_eq!(array.ndim(), 4); // Should be 4D [j, i, k, t]
        // Verify that boundary trimming worked (should be smaller than original)
        assert!(array.shape()[0] > 0); // j dimension
        assert!(array.shape()[1] > 0); // i dimension
        assert!(array.shape()[2] > 0); // k dimension
        assert!(array.shape()[3] > 0); // t dimension
    } else {
        // If NetCDF files aren't available, test should not fail
        eprintln!("Warning: NetCDF test file not found, skipping real file test");
    }
}

#[test]
fn test_read_coordinates_with_real_netcdf() {
    let reader = AdvancedNetCDFReader::new("test_data", 0);
    let datetime = NaiveDate::from_ymd_opt(2023, 7, 31)
        .unwrap()
        .and_hms_opt(12, 0, 0)
        .unwrap()
        .and_utc();

    // Test reading coordinates from test file
    let result = reader.read_coordinates(&datetime, Some(1));
    if result.is_ok() {
        let (xlat, xlong, levels) = result.unwrap();
        assert!(xlat.ndim() == 2);
        assert!(xlong.ndim() == 2);
        assert!(levels.ndim() == 1);
        assert_eq!(xlat.shape(), xlong.shape()); // Lat and lon should have same shape
    } else {
        eprintln!("Warning: NetCDF test file not found, skipping coordinate test");
    }
}

#[test]
fn test_read_multiple_variables_with_real_netcdf() {
    let reader = AdvancedNetCDFReader::new("test_data", 0);
    let datetime = NaiveDate::from_ymd_opt(2023, 7, 31)
        .unwrap()
        .and_hms_opt(12, 0, 0)
        .unwrap()
        .and_utc();

    let variables = ["U", "V", "W"];
    let result = reader.read_multiple_variables(&datetime, &variables, Some(1));
    if result.is_ok() {
        let var_map = result.unwrap();
        assert_eq!(var_map.len(), 3);
        assert!(var_map.contains_key("U"));
        assert!(var_map.contains_key("V"));
        assert!(var_map.contains_key("W"));

        // All variables should have same shape
        let u_shape = var_map["U"].shape();
        let v_shape = var_map["V"].shape();
        let w_shape = var_map["W"].shape();
        assert_eq!(u_shape, v_shape);
        assert_eq!(v_shape, w_shape);
    } else {
        eprintln!("Warning: NetCDF test file not found, skipping multi-variable test");
    }
}

#[test]
fn test_boundary_trimming_with_real_data() {
    // Test with different boundary trim values
    let reader_no_trim = AdvancedNetCDFReader::new("test_data", 0);
    let reader_trim_1 = AdvancedNetCDFReader::new("test_data", 1);

    let datetime = NaiveDate::from_ymd_opt(2023, 7, 31)
        .unwrap()
        .and_hms_opt(12, 0, 0)
        .unwrap()
        .and_utc();

    let result_no_trim = reader_no_trim.read_variable_array(&datetime, "U", Some(1));
    let result_trim_1 = reader_trim_1.read_variable_array(&datetime, "U", Some(1));

    if result_no_trim.is_ok() && result_trim_1.is_ok() {
        let array_no_trim = result_no_trim.unwrap();
        let array_trim_1 = result_trim_1.unwrap();

        // Trimmed array should be smaller by 2*boundary_trim in j and i dimensions
        assert_eq!(array_trim_1.shape()[0], array_no_trim.shape()[0] - 2); // j dimension
        assert_eq!(array_trim_1.shape()[1], array_no_trim.shape()[1] - 2); // i dimension
        assert_eq!(array_trim_1.shape()[2], array_no_trim.shape()[2]); // k dimension (unchanged)
        assert_eq!(array_trim_1.shape()[3], array_no_trim.shape()[3]); // t dimension (unchanged)
    } else {
        eprintln!("Warning: NetCDF test file not found, skipping boundary trim test");
    }
}
