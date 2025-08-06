use ndarray::{Array1, Array2, Array3, Array4, Zip};
use num_traits::Float;

/// Generic linear interpolation between two values
/// Replicates Fortran lin_interp with generic Float type
pub fn lin_interp<T: Float>(v0: T, v1: T, fac: T) -> T {
    v0 + (v1 - v0) * fac
}

/// Linear interpolation between two points (traditional interface)
pub fn linear_interpolate(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    if (x1 - x0).abs() < f64::EPSILON {
        return y0; // Avoid division by zero
    }
    let fac = (x - x0) / (x1 - x0);
    lin_interp(y0, y1, fac)
}

/// 3D linear interpolation along one axis
/// Replicates Fortran lin_interp3D functionality
pub fn lin_interp3d<T: Float + Copy>(v0: T, v1: T, fac: T) -> T {
    lin_interp(v0, v1, fac)
}

/// Parallel linear interpolation for ndarray
/// Uses Zip + par_azip! for SIMD/parallel processing
pub fn lin_interp_array<T>(v0: &Array1<T>, v1: &Array1<T>, fac: T) -> Array1<T>
where
    T: Float + Copy + Send + Sync,
{
    let mut result = Array1::zeros(v0.len());

    Zip::from(&mut result)
        .and(v0)
        .and(v1)
        .par_for_each(|res, &val0, &val1| {
            *res = lin_interp(val0, val1, fac);
        });

    result
}

/// Parallel linear interpolation for 2D arrays
pub fn lin_interp_array2d<T>(v0: &Array2<T>, v1: &Array2<T>, fac: T) -> Array2<T>
where
    T: Float + Copy + Send + Sync,
{
    let mut result = Array2::zeros(v0.dim());

    Zip::from(&mut result)
        .and(v0)
        .and(v1)
        .par_for_each(|res, &val0, &val1| {
            *res = lin_interp(val0, val1, fac);
        });

    result
}

/// Parallel linear interpolation for 3D arrays
pub fn lin_interp_array3d<T>(v0: &Array3<T>, v1: &Array3<T>, fac: T) -> Array3<T>
where
    T: Float + Copy + Send + Sync,
{
    let mut result = Array3::zeros(v0.dim());

    Zip::from(&mut result)
        .and(v0)
        .and(v1)
        .par_for_each(|res, &val0, &val1| {
            *res = lin_interp(val0, val1, fac);
        });

    result
}

/// Generic bilinear interpolation (replicates Fortran bilin_interp)
/// Interpolates between four corner values using two interpolation factors
pub fn bilin_interp<T: Float>(f00: T, f01: T, f10: T, f11: T, fac_x: T, fac_y: T) -> T {
    let f_y0 = lin_interp(f00, f10, fac_x);
    let f_y1 = lin_interp(f01, f11, fac_x);
    lin_interp(f_y0, f_y1, fac_y)
}

/// Bilinear interpolation in 2D grid (traditional interface)
#[allow(clippy::too_many_arguments)]
pub fn bilinear_interpolate(
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64, // Grid coordinates
    f00: f64,
    f01: f64,
    f10: f64,
    f11: f64, // Function values at grid corners
    x: f64,
    y: f64, // Interpolation point
) -> f64 {
    let dx = x1 - x0;
    let dy = y1 - y0;

    if dx.abs() < f64::EPSILON || dy.abs() < f64::EPSILON {
        return f00; // Degenerate case
    }

    let fac_x = (x - x0) / dx;
    let fac_y = (y - y0) / dy;

    bilin_interp(f00, f01, f10, f11, fac_x, fac_y)
}

/// Trilinear interpolation in 3D grid
#[allow(clippy::too_many_arguments)]
pub fn trilinear_interpolate(
    x0: f64,
    y0: f64,
    z0: f64,
    x1: f64,
    y1: f64,
    z1: f64, // Grid coordinates
    f000: f64,
    f001: f64,
    f010: f64,
    f011: f64, // Function values at 8 corners
    f100: f64,
    f101: f64,
    f110: f64,
    f111: f64,
    x: f64,
    y: f64,
    z: f64, // Interpolation point
) -> f64 {
    let dx = x1 - x0;
    let dy = y1 - y0;
    let dz = z1 - z0;

    if dx.abs() < f64::EPSILON || dy.abs() < f64::EPSILON || dz.abs() < f64::EPSILON {
        return f000; // Degenerate case
    }

    let _wx = (x - x0) / dx;
    let _wy = (y - y0) / dy;
    let wz = (z - z0) / dz;

    // Interpolate along z-axis first
    let f00 = f000 * (1.0 - wz) + f001 * wz;
    let f01 = f010 * (1.0 - wz) + f011 * wz;
    let f10 = f100 * (1.0 - wz) + f101 * wz;
    let f11 = f110 * (1.0 - wz) + f111 * wz;

    // Then bilinear interpolation in xy-plane
    bilinear_interpolate(x0, y0, x1, y1, f00, f01, f10, f11, x, y)
}

/// 4D interpolation (space + time) for meteorological data
#[allow(clippy::too_many_arguments)]
pub fn quadrilinear_interpolate(
    lon: f64,
    lat: f64,
    lev: f64,
    time: f64, // Interpolation point
    lon0: f64,
    lon1: f64, // Longitude bounds
    lat0: f64,
    lat1: f64, // Latitude bounds
    lev0: f64,
    lev1: f64, // Level bounds
    time0: f64,
    time1: f64,                // Time bounds
    data: &[[[f64; 2]; 2]; 2], // Data array [time][level][lat][lon]
) -> f64 {
    // Get weights
    let _wlon = if (lon1 - lon0).abs() < f64::EPSILON {
        0.0
    } else {
        (lon - lon0) / (lon1 - lon0)
    };
    let _wlat = if (lat1 - lat0).abs() < f64::EPSILON {
        0.0
    } else {
        (lat - lat0) / (lat1 - lat0)
    };
    let wlev = if (lev1 - lev0).abs() < f64::EPSILON {
        0.0
    } else {
        (lev - lev0) / (lev1 - lev0)
    };
    let wtime = if (time1 - time0).abs() < f64::EPSILON {
        0.0
    } else {
        (time - time0) / (time1 - time0)
    };

    // Interpolate in space for each time level
    let val_t0_l0 = bilinear_interpolate(
        lon0,
        lat0,
        lon1,
        lat1,
        data[0][0][0],
        data[0][0][1],
        data[0][1][0],
        data[0][1][1],
        lon,
        lat,
    );
    let val_t0_l1 = bilinear_interpolate(
        lon0,
        lat0,
        lon1,
        lat1,
        data[0][1][0],
        data[0][1][1],
        data[0][1][0],
        data[0][1][1],
        lon,
        lat,
    );
    let val_t1_l0 = bilinear_interpolate(
        lon0,
        lat0,
        lon1,
        lat1,
        data[1][0][0],
        data[1][0][1],
        data[1][1][0],
        data[1][1][1],
        lon,
        lat,
    );
    let val_t1_l1 = bilinear_interpolate(
        lon0,
        lat0,
        lon1,
        lat1,
        data[1][1][0],
        data[1][1][1],
        data[1][1][0],
        data[1][1][1],
        lon,
        lat,
    );

    // Interpolate in vertical
    let val_t0 = val_t0_l0 * (1.0 - wlev) + val_t0_l1 * wlev;
    let val_t1 = val_t1_l0 * (1.0 - wlev) + val_t1_l1 * wlev;

    // Interpolate in time
    val_t0 * (1.0 - wtime) + val_t1 * wtime
}

/// Find grid indices and weights for interpolation
pub fn find_grid_indices(coords: &[f64], target: f64) -> Result<(usize, usize, f64), String> {
    if coords.is_empty() {
        return Err("Empty coordinate array".to_string());
    }

    // Handle extrapolation cases
    if target <= coords[0] {
        return Ok((0, 0, 0.0));
    }
    if target >= coords[coords.len() - 1] {
        let last = coords.len() - 1;
        return Ok((last, last, 0.0));
    }

    // Binary search for insertion point
    let mut left = 0;
    let mut right = coords.len() - 1;

    while right - left > 1 {
        let mid = (left + right) / 2;
        if coords[mid] <= target {
            left = mid;
        } else {
            right = mid;
        }
    }

    // Calculate weight
    let weight = (target - coords[left]) / (coords[right] - coords[left]);

    Ok((left, right, weight))
}

/// Type alias for complex coordinate results
type CoordinatesResult = Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), String>;

/// Generic trait for data access that works with different readers (NetCDF, Zarr, etc.)
/// This trait abstracts away the specific data format and provides a common interface
/// for meteorological field interpolation.
pub trait FieldDataAccess {
    /// Get the shape of the data field [time, level, lat, lon]
    fn get_shape(&self) -> (usize, usize, usize, usize);
    
    /// Get a data value at specific indices [time, level, lat, lon]
    fn get_value(&self, time_idx: usize, level_idx: usize, lat_idx: usize, lon_idx: usize) -> Result<f64, String>;
    
    /// Get coordinate arrays
    fn get_coordinates(&self) -> CoordinatesResult; // (lon, lat, lev, time)
}

/// Implementation for legacy nested Vec structure (for backward compatibility)
pub struct VecFieldData<'a> {
    pub data: &'a [Vec<Vec<Vec<f64>>>],
    pub longitudes: &'a [f64],
    pub latitudes: &'a [f64], 
    pub levels: &'a [f64],
    pub times: &'a [f64],
}

impl FieldDataAccess for VecFieldData<'_> {
    fn get_shape(&self) -> (usize, usize, usize, usize) {
        let nt = self.data.len();
        let nk = if nt > 0 { self.data[0].len() } else { 0 };
        let nj = if nk > 0 { self.data[0][0].len() } else { 0 };
        let ni = if nj > 0 { self.data[0][0][0].len() } else { 0 };
        (nt, nk, nj, ni)
    }
    
    fn get_value(&self, time_idx: usize, level_idx: usize, lat_idx: usize, lon_idx: usize) -> Result<f64, String> {
        if time_idx >= self.data.len() {
            return Err(format!("Time index {} out of bounds", time_idx));
        }
        if level_idx >= self.data[time_idx].len() {
            return Err(format!("Level index {} out of bounds", level_idx));
        }
        if lat_idx >= self.data[time_idx][level_idx].len() {
            return Err(format!("Latitude index {} out of bounds", lat_idx));
        }
        if lon_idx >= self.data[time_idx][level_idx][lat_idx].len() {
            return Err(format!("Longitude index {} out of bounds", lon_idx));
        }
        Ok(self.data[time_idx][level_idx][lat_idx][lon_idx])
    }
    
    fn get_coordinates(&self) -> CoordinatesResult {
        Ok((self.longitudes.to_vec(), self.latitudes.to_vec(), self.levels.to_vec(), self.times.to_vec()))
    }
}

/// Implementation for ndarray-based data (modern approach)
pub struct ArrayFieldData<'a> {
    pub data: &'a Array4<f32>, // [j, i, k, t] layout from readers
    pub longitudes: &'a [f64],
    pub latitudes: &'a [f64],
    pub levels: &'a [f64], 
    pub times: &'a [f64],
}

impl FieldDataAccess for ArrayFieldData<'_> {
    fn get_shape(&self) -> (usize, usize, usize, usize) {
        let shape = self.data.shape();
        // Data is in [j, i, k, t] layout, convert to [t, k, j, i] for API compatibility
        (shape[3], shape[2], shape[0], shape[1])
    }
    
    fn get_value(&self, time_idx: usize, level_idx: usize, lat_idx: usize, lon_idx: usize) -> Result<f64, String> {
        let shape = self.data.shape();
        
        // Convert from [t, k, j, i] indices to [j, i, k, t] indices
        if lat_idx >= shape[0] || lon_idx >= shape[1] || level_idx >= shape[2] || time_idx >= shape[3] {
            return Err(format!("Index out of bounds: t={}, k={}, j={}, i={}", time_idx, level_idx, lat_idx, lon_idx));
        }
        
        // Access data with [j, i, k, t] indexing
        let value = self.data[[lat_idx, lon_idx, level_idx, time_idx]];
        Ok(value as f64)
    }
    
    fn get_coordinates(&self) -> CoordinatesResult {
        Ok((self.longitudes.to_vec(), self.latitudes.to_vec(), self.levels.to_vec(), self.times.to_vec()))
    }
}

/// Generic interpolation function that works with any FieldDataAccess implementation
pub fn interpolate_field_generic<T: FieldDataAccess>(
    field_accessor: &T,
    target_lon: f64,
    target_lat: f64,
    target_lev: f64,
    target_time: f64,
) -> Result<f64, String> {
    let (longitudes, latitudes, levels, times) = field_accessor.get_coordinates()?;
    
    // Find indices and weights for each dimension
    let (lon_i0, lon_i1, lon_w) = find_grid_indices(&longitudes, target_lon)?;
    let (lat_i0, lat_i1, lat_w) = find_grid_indices(&latitudes, target_lat)?;
    let (lev_i0, lev_i1, lev_w) = find_grid_indices(&levels, target_lev)?;
    let (time_i0, time_i1, time_w) = find_grid_indices(&times, target_time)?;

    // Extract 2x2x2x2 data cube
    let mut data_cube = [[[0.0; 2]; 2]; 2];

    for (t, data_t) in data_cube.iter_mut().enumerate() {
        let ti = if t == 0 { time_i0 } else { time_i1 };
        for (l, data_tl) in data_t.iter_mut().enumerate() {
            let li = if l == 0 { lev_i0 } else { lev_i1 };
            for (lat, data_value) in data_tl.iter_mut().enumerate() {
                let lati = if lat == 0 { lat_i0 } else { lat_i1 };
                let val0 = field_accessor.get_value(ti, li, lati, lon_i0)?;
                let val1 = field_accessor.get_value(ti, li, lati, lon_i1)?;
                *data_value = val0 * (1.0 - lon_w) + val1 * lon_w;
            }
        }
    }

    // Interpolate in remaining dimensions
    let mut result = 0.0;
    for (t, data_t) in data_cube.iter().enumerate() {
        let tw = if t == 0 { 1.0 - time_w } else { time_w };
        for (l, data_tl) in data_t.iter().enumerate() {
            let lw = if l == 0 { 1.0 - lev_w } else { lev_w };
            let val0 = data_tl[0];
            let val1 = data_tl[1];
            let val_interp = val0 * (1.0 - lat_w) + val1 * lat_w;
            result += val_interp * tw * lw;
        }
    }

    Ok(result)
}

/// Legacy interpolation function for backward compatibility
#[allow(clippy::too_many_arguments)]
pub fn interpolate_meteo_field(
    field_data: &[Vec<Vec<Vec<f64>>>], // [time][level][lat][lon]
    longitudes: &[f64],
    latitudes: &[f64],
    levels: &[f64],
    times: &[f64],
    target_lon: f64,
    target_lat: f64,
    target_lev: f64,
    target_time: f64,
) -> Result<f64, String> {
    // Use the generic interpolation with legacy Vec wrapper
    let field_accessor = VecFieldData {
        data: field_data,
        longitudes,
        latitudes, 
        levels,
        times,
    };
    
    interpolate_field_generic(&field_accessor, target_lon, target_lat, target_lev, target_time)
}
