/// Common functionality shared across data loaders
use super::ReaderError;
use chrono::{DateTime, Utc};
use ndarray::Array4;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Cache key generation for meteorological data
pub fn generate_cache_key(datetime: &DateTime<Utc>, prefix: Option<&str>) -> String {
    match prefix {
        Some(p) => format!("{}_{}", p, datetime.format("%Y-%m-%d_%H:%M:%S")),
        None => datetime.format("%Y-%m-%d_%H:%M:%S").to_string(),
    }
}

/// Apply boundary trimming to a 4D array with layout [j, i, k, t]
pub fn apply_boundary_trimming(array: Array4<f32>, boundary_trim: usize) -> Array4<f32> {
    if boundary_trim > 0 {
        let bdy = boundary_trim;
        let shape = array.shape();
        let (nj, ni) = (shape[0], shape[1]);

        array
            .slice(ndarray::s![bdy..(nj - bdy), bdy..(ni - bdy), .., ..])
            .to_owned()
    } else {
        array
    }
}

/// Build NetCDF filename following WRF conventions
/// This ports the get_filename logic typically used in meteorological models
pub fn get_filename(
    base_path: &str,
    datetime: &DateTime<Utc>,
    domain: Option<u32>,
    file_type: &str,
) -> PathBuf {
    let domain_suffix = match domain {
        Some(d) => format!("_d{:02}", d),
        None => String::new(),
    };

    let filename = format!(
        "{}{}_{}.nc",
        file_type,
        domain_suffix,
        datetime.format("%Y-%m-%d_%H:%M:%S")
    );

    Path::new(base_path).join(filename)
}

/// Open NetCDF file and validate it exists
pub fn open_netcdf_file(
    base_path: &str,
    datetime: &DateTime<Utc>,
    domain: Option<u32>,
    file_type: &str,
) -> Result<netcdf::File, ReaderError> {
    let filename = get_filename(base_path, datetime, domain, file_type);

    if !filename.exists() {
        return Err(ReaderError::FileNotFound(
            filename.to_string_lossy().to_string(),
        ));
    }

    Ok(netcdf::open(&filename)?)
}

/// Generic cache retrieval function
pub fn get_from_cache<T: Clone>(
    cache: &HashMap<String, T>,
    cache_key: &str,
) -> Option<T> {
    cache.get(cache_key).cloned()
}

/// Generic cache insertion function
pub fn insert_into_cache<T>(
    cache: &mut HashMap<String, T>,
    cache_key: String,
    data: T,
) {
    cache.insert(cache_key, data);
}

/// Handle missing variables in data loading with warning and optional fallback
pub fn handle_missing_variable(
    var_name: &str,
    error: &ReaderError,
    fallback_data: Option<&HashMap<String, Array4<f32>>>,
) -> Option<Array4<f32>> {
    eprintln!("Warning: Could not load variable {}: {}", var_name, error);
    
    // Create zeros array if we have reference data
    if let Some(data) = fallback_data {
        if let Some(ref_array) = data.values().next() {
            return Some(Array4::zeros(ref_array.raw_dim()));
        }
    }
    
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;
    use ndarray::Array4;

    #[test]
    fn test_generate_cache_key() {
        let dt = Utc.with_ymd_and_hms(2023, 1, 15, 12, 0, 0).unwrap();
        
        let key = generate_cache_key(&dt, None);
        assert_eq!(key, "2023-01-15_12:00:00");
        
        let prefixed_key = generate_cache_key(&dt, Some("test"));
        assert_eq!(prefixed_key, "test_2023-01-15_12:00:00");
    }
    
    #[test]
    fn test_apply_boundary_trimming() {
        let array = Array4::<f32>::zeros((10, 10, 5, 1));
        
        // No trimming
        let result = apply_boundary_trimming(array.clone(), 0);
        assert_eq!(result.shape(), &[10, 10, 5, 1]);
        
        // With trimming
        let result = apply_boundary_trimming(array, 2);
        assert_eq!(result.shape(), &[6, 6, 5, 1]);
    }
    
    #[test]
    fn test_get_filename() {
        let dt = Utc.with_ymd_and_hms(2023, 1, 15, 12, 0, 0).unwrap();
        
        let filename = get_filename("/data", &dt, Some(1), "wrfout");
        assert_eq!(
            filename.to_string_lossy(),
            "/data/wrfout_d01_2023-01-15_12:00:00.nc"
        );
        
        let filename_no_domain = get_filename("/data", &dt, None, "wrfout");
        assert_eq!(
            filename_no_domain.to_string_lossy(),
            "/data/wrfout_2023-01-15_12:00:00.nc"
        );
    }
}
