/// Utility functions for NetCDF data processing
use std::collections::HashMap;

/// Standard meteorological variable names used in WRF NetCDF files
pub const STANDARD_METEO_VARS: &[&str] = &["U", "V", "W", "T", "QVAPOR", "RAIN"];

/// Get the appropriate SI units for a meteorological variable
pub fn get_variable_units(variable_name: &str) -> String {
    match variable_name {
        // Wind components
        "U" | "V" | "W" | "u_wind" | "v_wind" | "w_wind" => "m/s".to_string(),
        // Temperature
        "T" | "temperature" => "K".to_string(),
        // Water vapor mixing ratio
        "QVAPOR" => "kg/kg".to_string(),
        // Precipitation
        "RAIN" | "RAINC" | "RAINNC" => "mm".to_string(),
        // Pressure
        "P" | "PB" | "pressure" => "Pa".to_string(),
        // Geopotential
        "PH" | "PHB" | "geopotential" => "m²/s²".to_string(),
        // Density
        "RHO" | "density" => "kg/m³".to_string(),
        // Default for unknown variables
        _ => "unknown".to_string(),
    }
}

/// Get the missing value commonly used for a meteorological variable
pub fn get_missing_value(variable_name: &str) -> f32 {
    match variable_name {
        // Standard WRF missing value
        _ => -9999.0,
    }
}

/// Create a map of variable names to their units
pub fn create_units_map() -> HashMap<String, String> {
    let mut units = HashMap::new();

    for &var in STANDARD_METEO_VARS {
        units.insert(var.to_string(), get_variable_units(var));
    }

    // Add coordinate variables
    units.insert("XLAT".to_string(), "degrees_north".to_string());
    units.insert("XLONG".to_string(), "degrees_east".to_string());
    units.insert(
        "XTIME".to_string(),
        "minutes since simulation start".to_string(),
    );

    units
}

/// Check if a variable name represents a standard meteorological field
pub fn is_standard_meteo_var(variable_name: &str) -> bool {
    STANDARD_METEO_VARS.contains(&variable_name)
}

/// Check if a variable name represents a coordinate variable
pub fn is_coordinate_var(variable_name: &str) -> bool {
    matches!(variable_name, "XLAT" | "XLONG" | "XTIME" | "Times")
}

/// Validate NetCDF dimension ordering for WRF files
/// Expected order: [Time, bottom_top, south_north, west_east]
pub fn validate_wrf_dimensions(shape: &[usize]) -> Result<(), String> {
    if shape.len() != 4 {
        return Err(format!("Expected 4 dimensions, got {}", shape.len()));
    }

    // Basic sanity checks
    let (nt, nk, nj, ni) = (shape[0], shape[1], shape[2], shape[3]);

    if nt == 0 || nk == 0 || nj == 0 || ni == 0 {
        return Err("All dimensions must be greater than 0".to_string());
    }

    if nj < 10 || ni < 10 {
        return Err("Spatial dimensions seem too small for a realistic grid".to_string());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_variable_units() {
        assert_eq!(get_variable_units("U"), "m/s");
        assert_eq!(get_variable_units("T"), "K");
        assert_eq!(get_variable_units("QVAPOR"), "kg/kg");
        assert_eq!(get_variable_units("unknown_var"), "unknown");
    }

    #[test]
    fn test_is_standard_meteo_var() {
        assert!(is_standard_meteo_var("U"));
        assert!(is_standard_meteo_var("T"));
        assert!(!is_standard_meteo_var("XLAT"));
    }

    #[test]
    fn test_is_coordinate_var() {
        assert!(is_coordinate_var("XLAT"));
        assert!(is_coordinate_var("XLONG"));
        assert!(!is_coordinate_var("U"));
    }

    #[test]
    fn test_validate_wrf_dimensions() {
        assert!(validate_wrf_dimensions(&[3, 5, 20, 30]).is_ok());
        assert!(validate_wrf_dimensions(&[3, 5, 20]).is_err()); // Too few dimensions
        assert!(validate_wrf_dimensions(&[3, 5, 5, 30]).is_err()); // Too small spatial grid
    }
}
