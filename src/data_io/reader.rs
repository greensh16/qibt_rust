// Removed utils import
use super::{MeteoData, MeteoField};
use chrono::{DateTime, Utc, Duration};
use ndarray::{Array1, Array2, Array4};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use thiserror::Error;
use crate::config::Config;
use crate::io::{DataReader, DataReaderError, VariableInfo, DimensionInfo, Attributes, AttributeValue, CoordinateArrays};

/// Type alias for coordinate result to reduce complexity
type CoordinateResult = (Array2<f32>, Array2<f32>, Array1<f32>);

/// Common interface for data loaders to guarantee identical external behavior
pub trait DataLoader {
    /// Load standard meteorological data for a given time
    /// This matches the get_data() function from Fortran QIBT
    fn get_data(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<HashMap<String, Array4<f32>>, ReaderError>;

    /// Load cloud/precipitation mixing ratios
    /// This matches the get_data_mixtot() function from Fortran QIBT
    fn get_data_mixtot(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<HashMap<String, Array4<f32>>, ReaderError>;

    /// Load precipitation and latent heat data
    fn get_precipitation_data(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<HashMap<String, Array4<f32>>, ReaderError>;

    /// Load surface pressure data
    fn get_surface_pressure(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<Array4<f32>, ReaderError>;

    /// Load data for entire simulation period
    /// This provides the multi-day, multi-file loading capability from Fortran QIBT
    fn load_simulation_data(
        &mut self,
        config: &Config,
    ) -> Result<SimulationDataset, ReaderError>;

    /// Clear the data cache to free memory
    fn clear_cache(&mut self);

    /// Generate WRF filename for specific datetime and file type
    /// Matches the get_filename() function from Fortran QIBT
    fn get_wrfout_filename(&self, datetime: &DateTime<Utc>) -> PathBuf;

    /// Get a slice of variable data for interpolation purposes
    /// This method enables efficient access to subsets of meteorological data
    fn get_variable_slice(
        &mut self,
        datetime: &DateTime<Utc>,
        variable_name: &str,
        slice_params: Option<(usize, usize, usize, usize)>, // (j_start, j_end, i_start, i_end)
    ) -> Result<Array4<f32>, ReaderError>;

    /// Cache data for a specific time step to improve performance
    /// This method allows pre-loading data for anticipated access patterns
    fn cache_step(&mut self, datetime: &DateTime<Utc>) -> Result<(), ReaderError>;
}

/// NetCDF reader for meteorological data
pub struct NetCDFReader {
    pub file_path: String,
}

#[derive(Error, Debug)]
pub enum ReaderError {
    #[error("NetCDF error: {0}")]
    Netcdf(#[from] netcdf::Error),

    #[error("Variable not found: {0}")]
    MissingVariable(String),

    #[error("Data conversion error")]
    ConversionError,

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid datetime format: {0}")]
    DateTimeError(String),

    #[error("File not found: {0}")]
    FileNotFound(String),
}

impl NetCDFReader {
    /// Create a new NetCDF reader
    pub fn new(file_path: impl AsRef<Path>) -> Self {
        Self {
            file_path: file_path.as_ref().to_string_lossy().to_string(),
        }
    }

    /// Read complete meteorological dataset
    pub fn read_meteo_data(&self) -> Result<MeteoData, ReaderError> {
        // For now, return a placeholder implementation
        // This would be implemented when NetCDF files are available
        let mut data = MeteoData::new();

        // Read the grid information (placeholder)
        for var_name in &["U", "V", "W", "T", "QVAPOR", "RAIN"] {
            // Add to the data fields
            data.fields.push(MeteoField {
                name: var_name.to_string(),
                data: vec![], // Placeholder: Convert Array4 to Vec
                missing_value: -9999.0,
                units: match *var_name {
                    "U" | "V" | "W" => "m/s".to_string(),
                    "T" => "K".to_string(),
                    "QVAPOR" => "kg/kg".to_string(),
                    "RAIN" => "mm".to_string(),
                    _ => "unknown".to_string(),
                },
            });
        }

        println!("Reading meteorological data from: {}", self.file_path);
        Ok(data)
    }

    /// Read specific meteorological field
    pub fn read_field(&self, field_name: &str) -> Result<MeteoField, String> {
        // TODO: Implement field-specific reading
        let field = MeteoField {
            name: field_name.to_string(),
            data: vec![], // Placeholder
            missing_value: -9999.0,
            units: match field_name {
                "temperature" => "K".to_string(),
                "u_wind" | "v_wind" | "w_wind" => "m/s".to_string(),
                "pressure" => "Pa".to_string(),
                "geopotential" => "m²/s²".to_string(),
                _ => "unknown".to_string(),
            },
        };

        println!("Reading field: {}", field_name);
        Ok(field)
    }

    /// Check if file exists and is readable
    pub fn validate_file(&self) -> Result<(), String> {
        let path = Path::new(&self.file_path);
        if !path.exists() {
            return Err(format!("File does not exist: {}", self.file_path));
        }
        if !path.is_file() {
            return Err(format!("Path is not a file: {}", self.file_path));
        }
        Ok(())
    }

    /// Get file metadata
    pub fn get_metadata(&self) -> Result<NetCDFFileMetadata, String> {
        // TODO: Implement metadata reading
        Ok(NetCDFFileMetadata {
            dimensions: vec![
                "time".to_string(),
                "level".to_string(),
                "lat".to_string(),
                "lon".to_string(),
            ],
            variables: vec![
                "temperature".to_string(),
                "u_wind".to_string(),
                "v_wind".to_string(),
            ],
            global_attributes: std::collections::HashMap::new(),
        })
    }
}

/// NetCDF-specific file metadata structure
#[derive(Debug)]
pub struct NetCDFFileMetadata {
    pub dimensions: Vec<String>,
    pub variables: Vec<String>,
    pub global_attributes: std::collections::HashMap<String, String>,
}

/// Read wind components at specific location and time
pub fn read_winds_at_point(
    _reader: &NetCDFReader,
    _lon: f64,
    _lat: f64,
    _pressure: f64,
    _time: f64,
) -> Result<(f64, f64, f64), String> {
    // TODO: Implement point interpolation from NetCDF data
    // This would involve reading the full fields and interpolating

    // Placeholder implementation
    let u_wind = 10.0; // m/s
    let v_wind = 5.0; // m/s
    let w_wind = 0.1; // m/s (vertical velocity)

    Ok((u_wind, v_wind, w_wind))
}

/// Read temperature at specific location and time
pub fn read_temperature_at_point(
    _reader: &NetCDFReader,
    _lon: f64,
    _lat: f64,
    _pressure: f64,
    _time: f64,
) -> Result<f64, String> {
    // TODO: Implement point interpolation
    let temperature = 273.15 + 20.0; // K (20°C)
    Ok(temperature)
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

/// Advanced NetCDF reader with full functionality
pub struct AdvancedNetCDFReader {
    pub base_path: String,
    pub boundary_trim: usize,
}

impl AdvancedNetCDFReader {
    /// Create new advanced reader
    pub fn new(base_path: impl AsRef<str>, boundary_trim: usize) -> Self {
        Self {
            base_path: base_path.as_ref().to_string(),
            boundary_trim,
        }
    }

    /// Apply boundary trimming to a 4D array with layout [j, i, k, t]
    fn apply_boundary_trimming(&self, array: Array4<f32>) -> Array4<f32> {
        apply_boundary_trimming(array, self.boundary_trim)
    }

    /// Open NetCDF file and validate it exists
    fn open_netcdf_file(
        &self,
        datetime: &DateTime<Utc>,
        domain: Option<u32>,
        file_type: &str,
    ) -> Result<netcdf::File, ReaderError> {
        open_netcdf_file(&self.base_path, datetime, domain, file_type)
    }

    /// Read a specific variable as Array4<f32> with proper dimension ordering [j, i, k, t]
    pub fn read_variable_array(
        &self,
        datetime: &DateTime<Utc>,
        variable_name: &str,
        domain: Option<u32>,
    ) -> Result<Array4<f32>, ReaderError> {
        // Open the NetCDF file
        let file = self.open_netcdf_file(datetime, domain, "wrfout")?;

        // Get the variable
        let var = file
            .variable(variable_name)
            .ok_or_else(|| ReaderError::MissingVariable(variable_name.to_string()))?;

        // Read the data - NetCDF dimensions are typically [Time, bottom_top, south_north, west_east]
        let raw_data: Vec<f32> = var.get_values(..)?;
        let shape = var.dimensions().iter().map(|d| d.len()).collect::<Vec<_>>();

        if shape.len() != 4 {
            return Err(ReaderError::ConversionError);
        }

        let (nt, nk, nj, ni) = (shape[0], shape[1], shape[2], shape[3]);

        // Create Array4 from raw data with original shape [t, k, j, i]
        let original_array = Array4::from_shape_vec((nt, nk, nj, ni), raw_data)
            .map_err(|_| ReaderError::ConversionError)?;

        // Transform to [j, i, k, t] for interpolation-ready layout
        let transformed = original_array.permuted_axes([2, 3, 1, 0]);

        // Apply boundary trimming and return
        Ok(self.apply_boundary_trimming(transformed))
    }

    /// Read multiple variables efficiently
    pub fn read_multiple_variables(
        &self,
        datetime: &DateTime<Utc>,
        variables: &[&str],
        domain: Option<u32>,
    ) -> Result<HashMap<String, Array4<f32>>, ReaderError> {
        let mut result = HashMap::new();

        for &var_name in variables {
            let array = self.read_variable_array(datetime, var_name, domain)?;
            result.insert(var_name.to_string(), array);
        }

        Ok(result)
    }

    /// Read standard meteorological variables (U, V, W, T, QVAPOR, RAIN)
    pub fn read_meteo_variables(
        &self,
        datetime: &DateTime<Utc>,
        domain: Option<u32>,
    ) -> Result<HashMap<String, Array4<f32>>, ReaderError> {
        let standard_vars = ["U", "V", "W", "T", "QVAPOR", "P", "PB", "PH", "PHB"];
        self.read_multiple_variables(datetime, &standard_vars, domain)
    }

    /// Read coordinate variables (XLAT, XLONG, etc.)
    pub fn read_coordinates(
        &self,
        datetime: &DateTime<Utc>,
        domain: Option<u32>,
    ) -> Result<CoordinateResult, ReaderError> {
        // Open the NetCDF file
        let file = self.open_netcdf_file(datetime, domain, "wrfout")?;

        // Read XLAT (latitude)
        let xlat_var = file
            .variable("XLAT")
            .ok_or_else(|| ReaderError::MissingVariable("XLAT".to_string()))?;
        let xlat_data: Vec<f32> = xlat_var.get_values(..)?;
        let xlat_shape = xlat_var
            .dimensions()
            .iter()
            .map(|d| d.len())
            .collect::<Vec<_>>();
        let xlat = Array2::from_shape_vec((xlat_shape[0], xlat_shape[1]), xlat_data)
            .map_err(|_| ReaderError::ConversionError)?;

        // Read XLONG (longitude)
        let xlong_var = file
            .variable("XLONG")
            .ok_or_else(|| ReaderError::MissingVariable("XLONG".to_string()))?;
        let xlong_data: Vec<f32> = xlong_var.get_values(..)?;
        let xlong_shape = xlong_var
            .dimensions()
            .iter()
            .map(|d| d.len())
            .collect::<Vec<_>>();
        let xlong = Array2::from_shape_vec((xlong_shape[0], xlong_shape[1]), xlong_data)
            .map_err(|_| ReaderError::ConversionError)?;

        // Read pressure levels from P variable if available
        let levels = if let Some(p_var) = file.variable("P") {
            let p_data: Vec<f32> = p_var.get_values(..)?;
            let p_shape = p_var
                .dimensions()
                .iter()
                .map(|d| d.len())
                .collect::<Vec<_>>();
            let p_array =
                Array4::from_shape_vec((p_shape[0], p_shape[1], p_shape[2], p_shape[3]), p_data)
                    .map_err(|_| ReaderError::ConversionError)?;
            // Extract pressure profile from first time step and first grid point
            let profile = p_array.slice(ndarray::s![0, .., 0, 0]).to_owned();
            profile.into_dimensionality::<ndarray::Ix1>().unwrap()
        } else {
            // Default pressure levels if not available
            Array1::from(vec![100000.0, 95000.0, 90000.0, 85000.0, 80000.0])
        };

        Ok((xlat, xlong, levels))
    }
}

/// Multi-file data loading system matching Fortran QIBT functionality
pub struct MultiFileDataLoader {
    pub base_path: String,
    pub boundary_trim: usize,
    pub domain: u32,
    /// Cache for loaded meteorological data
    pub data_cache: HashMap<String, HashMap<String, Array4<f32>>>,
}

impl MultiFileDataLoader {
    /// Create new multi-file data loader
    pub fn new(base_path: impl AsRef<str>, boundary_trim: usize, domain: u32) -> Self {
        Self {
            base_path: base_path.as_ref().to_string(),
            boundary_trim,
            domain,
            data_cache: HashMap::new(),
        }
    }

    /// Generate WRF filename for specific datetime and file type
    /// Matches the get_filename() function from Fortran QIBT
    pub fn get_wrfout_filename(&self, datetime: &DateTime<Utc>) -> PathBuf {
        get_filename(&self.base_path, datetime, Some(self.domain), "wrfout")
    }

    /// Load standard meteorological data for a given time
    /// This matches the get_data() function from Fortran QIBT
    pub fn get_data(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<HashMap<String, Array4<f32>>, ReaderError> {
        let cache_key = generate_cache_key(datetime, None);
        
        // Check cache first
        if let Some(cached_data) = get_from_cache(&self.data_cache, &cache_key) {
            return Ok(cached_data);
        }

        let reader = AdvancedNetCDFReader::new(&self.base_path, self.boundary_trim);
        
        // Load standard atmospheric variables: U, V, W, T, QVAPOR, P, PH
        let variables = [
            "U",      // U-component wind
            "V",      // V-component wind 
            "W",      // Vertical velocity
            "T",      // Perturbation potential temperature
            "QVAPOR", // Water vapor mixing ratio
            "P",      // Perturbation pressure
            "PB",     // Base state pressure
            "PH",     // Perturbation geopotential
            "PHB",    // Base state geopotential
        ];

        let mut data = HashMap::new();
        
        for &var_name in &variables {
            match reader.read_variable_array(datetime, var_name, Some(self.domain)) {
                Ok(array) => {
                    data.insert(var_name.to_string(), array);
                },
                Err(e) => {
                    eprintln!("Warning: Could not load variable {}: {}", var_name, e);
                    // Continue loading other variables even if one fails
                }
            }
        }

        // Cache the loaded data
        insert_into_cache(&mut self.data_cache, cache_key, data.clone());
        
        Ok(data)
    }

    /// Load cloud/precipitation mixing ratios
    /// This matches the get_data_mixtot() function from Fortran QIBT
    pub fn get_data_mixtot(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<HashMap<String, Array4<f32>>, ReaderError> {
        let cache_key = generate_cache_key(datetime, Some("mixtot"));
        
        // Check cache first
        if let Some(cached_data) = get_from_cache(&self.data_cache, &cache_key) {
            return Ok(cached_data);
        }

        let reader = AdvancedNetCDFReader::new(&self.base_path, self.boundary_trim);
        
        // Load cloud and precipitation mixing ratios
        let variables = [
            "QCLOUD",  // Cloud water mixing ratio
            "QRAIN",   // Rain water mixing ratio
            "QICE",    // Ice mixing ratio
            "QSNOW",   // Snow mixing ratio
            "QGRAUP",  // Graupel mixing ratio
        ];

        let mut data = HashMap::new();
        
        for &var_name in &variables {
            match reader.read_variable_array(datetime, var_name, Some(self.domain)) {
                Ok(array) => {
                    data.insert(var_name.to_string(), array);
                },
                Err(e) => {
                    eprintln!("Warning: Could not load mixing ratio {}: {}", var_name, e);
                    // For missing variables, create zeros array if we have other data
                    if let Some(ref_array) = data.values().next() {
                        let zeros = Array4::zeros(ref_array.raw_dim());
                        data.insert(var_name.to_string(), zeros);
                    }
                }
            }
        }

        // Cache the loaded data
        insert_into_cache(&mut self.data_cache, cache_key, data.clone());
        
        Ok(data)
    }

    /// Load precipitation and latent heat data
    pub fn get_precipitation_data(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<HashMap<String, Array4<f32>>, ReaderError> {
        let cache_key = generate_cache_key(datetime, Some("precip"));
        
        // Check cache first
        if let Some(cached_data) = get_from_cache(&self.data_cache, &cache_key) {
            return Ok(cached_data);
        }

        let reader = AdvancedNetCDFReader::new(&self.base_path, self.boundary_trim);
        
        // Load precipitation and latent heat variables
        let variables = [
            "RAINC",    // Accumulated convective precipitation
            "RAINNC",   // Accumulated non-convective precipitation
            "LATENTC",  // Convective latent heating
            "LATENTNC", // Non-convective latent heating
        ];

        let mut data = HashMap::new();
        
        for &var_name in &variables {
            match reader.read_variable_array(datetime, var_name, Some(self.domain)) {
                Ok(array) => {
                    data.insert(var_name.to_string(), array);
                },
                Err(e) => {
                    eprintln!("Warning: Could not load precipitation variable {}: {}", var_name, e);
                }
            }
        }

        // Cache the loaded data
        insert_into_cache(&mut self.data_cache, cache_key, data.clone());
        
        Ok(data)
    }

    /// Load surface pressure data
    pub fn get_surface_pressure(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<Array4<f32>, ReaderError> {
        let reader = AdvancedNetCDFReader::new(&self.base_path, self.boundary_trim);
        reader.read_variable_array(datetime, "PSFC", Some(self.domain))
    }

    /// Generate sequence of datetime objects for multi-day simulation
    /// This matches the time stepping logic from Fortran QIBT
    pub fn generate_time_sequence(
        start_date: DateTime<Utc>,
        end_date: DateTime<Utc>,
        time_step_hours: f64,
    ) -> Vec<DateTime<Utc>> {
        let mut times = Vec::new();
        let step_duration = Duration::seconds((time_step_hours * 3600.0) as i64);
        
        let mut current_time = start_date;
        while current_time <= end_date {
            times.push(current_time);
            current_time += step_duration;
        }
        
        times
    }

    /// Load data for entire simulation period
    /// This provides the multi-day, multi-file loading capability from Fortran QIBT
    pub fn load_simulation_data(
        &mut self,
        config: &Config,
    ) -> Result<SimulationDataset, ReaderError> {
        let start_datetime = DateTime::from_timestamp(config.start_time as i64, 0)
            .ok_or_else(|| ReaderError::DateTimeError("Invalid start time".to_string()))?;
        let end_datetime = DateTime::from_timestamp(config.end_time as i64, 0)
            .ok_or_else(|| ReaderError::DateTimeError("Invalid end time".to_string()))?;
        
        let time_sequence = Self::generate_time_sequence(
            start_datetime,
            end_datetime,
            config.time_step / 3600.0, // Convert seconds to hours
        );

        let mut dataset = SimulationDataset::new();
        
        for datetime in &time_sequence {
            println!("Loading data for time: {}", datetime.format("%Y-%m-%d %H:%M:%S"));
            
            // Load main meteorological data
            let meteo_data = self.get_data(datetime)?;
            dataset.meteo_data.insert(*datetime, meteo_data);
            
            // Load mixing ratios
            if let Ok(mixtot_data) = self.get_data_mixtot(datetime) {
                dataset.mixtot_data.insert(*datetime, mixtot_data);
            }
            
            // Load precipitation data
            if let Ok(precip_data) = self.get_precipitation_data(datetime) {
                dataset.precip_data.insert(*datetime, precip_data);
            }
            
            // Load surface pressure
            if let Ok(surface_pressure) = self.get_surface_pressure(datetime) {
                dataset.surface_pressure.insert(*datetime, surface_pressure);
            }
        }
        
        // Load coordinates from the first time step
        if let Some(first_time) = time_sequence.first() {
            let reader = AdvancedNetCDFReader::new(&self.base_path, self.boundary_trim);
            if let Ok((xlat, xlong, levels)) = reader.read_coordinates(first_time, Some(self.domain)) {
                dataset.coordinates = Some((xlat, xlong, levels));
            }
        }
        
        println!("Loaded data for {} time steps", dataset.meteo_data.len());
        Ok(dataset)
    }

    /// Clear the data cache to free memory
    pub fn clear_cache(&mut self) {
        self.data_cache.clear();
    }

    /// Get a slice of variable data for interpolation purposes
    pub fn get_variable_slice(
        &mut self,
        datetime: &DateTime<Utc>,
        variable_name: &str,
        slice_params: Option<(usize, usize, usize, usize)>, // (j_start, j_end, i_start, i_end)
    ) -> Result<Array4<f32>, ReaderError> {
        let reader = AdvancedNetCDFReader::new(&self.base_path, self.boundary_trim);
        let full_array = reader.read_variable_array(datetime, variable_name, Some(self.domain))?;
        
        if let Some((j_start, j_end, i_start, i_end)) = slice_params {
            let sliced = full_array.slice(ndarray::s![j_start..j_end, i_start..i_end, .., ..]).to_owned();
            Ok(sliced)
        } else {
            Ok(full_array)
        }
    }

    /// Cache data for a specific time step to improve performance
    pub fn cache_step(&mut self, datetime: &DateTime<Utc>) -> Result<(), ReaderError> {
        // Pre-load all standard data for this time step
        let _ = self.get_data(datetime)?;
        let _ = self.get_data_mixtot(datetime);
        let _ = self.get_precipitation_data(datetime);
        let _ = self.get_surface_pressure(datetime);
        Ok(())
    }
}

/// Implementation of the DataLoader trait for MultiFileDataLoader
impl DataLoader for MultiFileDataLoader {
    fn get_data(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<HashMap<String, Array4<f32>>, ReaderError> {
        self.get_data(datetime)
    }

    fn get_data_mixtot(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<HashMap<String, Array4<f32>>, ReaderError> {
        self.get_data_mixtot(datetime)
    }

    fn get_precipitation_data(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<HashMap<String, Array4<f32>>, ReaderError> {
        self.get_precipitation_data(datetime)
    }

    fn get_surface_pressure(
        &mut self,
        datetime: &DateTime<Utc>,
    ) -> Result<Array4<f32>, ReaderError> {
        self.get_surface_pressure(datetime)
    }

    fn load_simulation_data(
        &mut self,
        config: &Config,
    ) -> Result<SimulationDataset, ReaderError> {
        self.load_simulation_data(config)
    }

    fn clear_cache(&mut self) {
        self.clear_cache()
    }

    fn get_wrfout_filename(&self, datetime: &DateTime<Utc>) -> PathBuf {
        self.get_wrfout_filename(datetime)
    }

    fn get_variable_slice(
        &mut self,
        datetime: &DateTime<Utc>,
        variable_name: &str,
        slice_params: Option<(usize, usize, usize, usize)>,
    ) -> Result<Array4<f32>, ReaderError> {
        self.get_variable_slice(datetime, variable_name, slice_params)
    }

    fn cache_step(&mut self, datetime: &DateTime<Utc>) -> Result<(), ReaderError> {
        self.cache_step(datetime)
    }
}

/// Complete dataset for simulation matching Fortran QIBT data structures
#[derive(Debug, Default)]
pub struct SimulationDataset {
    /// Meteorological data indexed by time
    pub meteo_data: HashMap<DateTime<Utc>, HashMap<String, Array4<f32>>>,
    /// Mixing ratio data indexed by time
    pub mixtot_data: HashMap<DateTime<Utc>, HashMap<String, Array4<f32>>>,
    /// Precipitation data indexed by time
    pub precip_data: HashMap<DateTime<Utc>, HashMap<String, Array4<f32>>>,
    /// Surface pressure indexed by time
    pub surface_pressure: HashMap<DateTime<Utc>, Array4<f32>>,
    /// Grid coordinates (lat, lon, pressure levels)
    pub coordinates: Option<(Array2<f32>, Array2<f32>, Array1<f32>)>,
}

impl SimulationDataset {
    pub fn new() -> Self {
        Self::default()
    }

    
    /// Get the available time steps
    pub fn get_time_steps(&self) -> Vec<DateTime<Utc>> {
        let mut times: Vec<_> = self.meteo_data.keys().copied().collect();
        times.sort();
        times
    }
    
    /// Get grid dimensions from the first available data
    pub fn get_grid_dimensions(&self) -> Option<(usize, usize, usize, usize)> {
        if let Some(data) = self.meteo_data.values().next() {
            if let Some(array) = data.values().next() {
                let shape = array.shape();
                return Some((shape[0], shape[1], shape[2], shape[3])); // (nj, ni, nk, nt)
            }
        }
        None
    }
}

// Helper functions for caching and NetCDF operations
fn generate_cache_key(datetime: &DateTime<Utc>, suffix: Option<&str>) -> String {
    match suffix {
        Some(s) => format!("{}-{}", datetime.format("%Y%m%d%H%M%S"), s),
        None => datetime.format("%Y%m%d%H%M%S").to_string(),
    }
}

fn get_from_cache(
    cache: &HashMap<String, HashMap<String, Array4<f32>>>,
    key: &str,
) -> Option<HashMap<String, Array4<f32>>> {
    cache.get(key).cloned()
}

fn insert_into_cache(
    cache: &mut HashMap<String, HashMap<String, Array4<f32>>>,
    key: String,
    data: HashMap<String, Array4<f32>>,
) {
    cache.insert(key, data);
}

fn apply_boundary_trimming(array: Array4<f32>, trim: usize) -> Array4<f32> {
    if trim == 0 {
        return array;
    }
    
    let shape = array.shape();
    if shape[0] <= 2 * trim || shape[1] <= 2 * trim {
        // If trimming would remove everything, return original
        return array;
    }
    
    // Trim boundaries: [j, i, k, t] -> [j-trim, i-trim, k, t]
    let trimmed = array.slice(ndarray::s![
        trim..shape[0]-trim,
        trim..shape[1]-trim,
        ..,
        ..
    ]).to_owned();
    
    trimmed
}

fn open_netcdf_file(
    base_path: &str,
    datetime: &DateTime<Utc>,
    domain: Option<u32>,
    file_type: &str,
) -> Result<netcdf::File, ReaderError> {
    let filename = get_filename(base_path, datetime, domain, file_type);
    
    if !filename.exists() {
        return Err(ReaderError::FileNotFound(format!(
            "NetCDF file not found: {}",
            filename.display()
        )));
    }
    
    netcdf::open(&filename).map_err(ReaderError::Netcdf)
}

/// Implementation of DataReader trait for NetCDFReader
impl DataReader for NetCDFReader {
    fn open(&mut self, path: &Path) -> Result<(), DataReaderError> {
        self.file_path = path.to_string_lossy().to_string();
        // Validate the file exists
        self.validate_file().map_err(DataReaderError::FileNotFound)?;
        Ok(())
    }
    
    fn is_open(&self) -> bool {
        Path::new(&self.file_path).exists()
    }
    
    fn close(&mut self) {
        // NetCDF files are stateless in this implementation
    }
    
    fn list_variables(&self) -> Result<Vec<String>, DataReaderError> {
        // Return standard meteorological variables
        Ok(vec![
            "U".to_string(),
            "V".to_string(), 
            "W".to_string(),
            "T".to_string(),
            "QVAPOR".to_string(),
            "P".to_string(),
            "PB".to_string(),
        ])
    }
    
    fn get_variable_info(&self, variable_name: &str) -> Result<VariableInfo, DataReaderError> {
        let info = VariableInfo {
            name: variable_name.to_string(),
            dimensions: vec!["time".to_string(), "level".to_string(), "lat".to_string(), "lon".to_string()],
            shape: vec![1, 10, 50, 50], // Placeholder dimensions
            dtype: "f32".to_string(),
            units: match variable_name {
                "U" | "V" | "W" => Some("m/s".to_string()),
                "T" => Some("K".to_string()),
                "QVAPOR" => Some("kg/kg".to_string()),
                "P" | "PB" => Some("Pa".to_string()),
                _ => None,
            },
            long_name: match variable_name {
                "U" => Some("U-component of wind".to_string()),
                "V" => Some("V-component of wind".to_string()),
                "W" => Some("Vertical velocity".to_string()),
                "T" => Some("Temperature".to_string()),
                "QVAPOR" => Some("Water vapor mixing ratio".to_string()),
                "P" => Some("Perturbation pressure".to_string()),
                "PB" => Some("Base state pressure".to_string()),
                _ => None,
            },
        };
        Ok(info)
    }
    
    fn list_dimensions(&self) -> Result<Vec<String>, DataReaderError> {
        Ok(vec![
            "time".to_string(),
            "level".to_string(),
            "lat".to_string(),
            "lon".to_string(),
        ])
    }
    
    fn get_dimension_info(&self, dimension_name: &str) -> Result<DimensionInfo, DataReaderError> {
        let info = DimensionInfo {
            name: dimension_name.to_string(),
            size: match dimension_name {
                "time" => 1,
                "level" => 10,
                "lat" => 50,
                "lon" => 50,
                _ => return Err(DataReaderError::MissingDimension(dimension_name.to_string())),
            },
            is_unlimited: dimension_name == "time",
        };
        Ok(info)
    }
    
    fn get_metadata(&self) -> Result<crate::io::FileMetadata, DataReaderError> {
        let dimensions = vec![
            DimensionInfo { name: "time".to_string(), size: 1, is_unlimited: true },
            DimensionInfo { name: "level".to_string(), size: 10, is_unlimited: false },
            DimensionInfo { name: "lat".to_string(), size: 50, is_unlimited: false },
            DimensionInfo { name: "lon".to_string(), size: 50, is_unlimited: false },
        ];
        
        let variables = self.list_variables()?
            .into_iter()
            .map(|name| self.get_variable_info(&name))
            .collect::<Result<Vec<_>, _>>()?;
        
        let global_attributes = Attributes::new();
        
        Ok(crate::io::FileMetadata {
            dimensions,
            variables,
            global_attributes,
        })
    }
    
    fn read_variable(&self, variable_name: &str) -> Result<Array4<f32>, DataReaderError> {
        // Create placeholder data with appropriate shape [j, i, k, t]
        let shape = (50, 50, 10, 1); // (lat, lon, level, time)
        let data = Array4::<f32>::zeros(shape);
        
        println!("NetCDF reader: reading variable {} with shape {:?}", variable_name, shape);
        Ok(data)
    }
    
    fn read_variable_slice(
        &self,
        variable_name: &str,
        _indices: &[(usize, usize, usize)]
    ) -> Result<Array4<f32>, DataReaderError> {
        // For simplicity, just return the full variable
        self.read_variable(variable_name)
    }
    
    fn read_variables(&self, variable_names: &[&str]) -> Result<HashMap<String, Array4<f32>>, DataReaderError> {
        let mut result = HashMap::new();
        
        for &name in variable_names {
            let data = self.read_variable(name)?;
            result.insert(name.to_string(), data);
        }
        
        Ok(result)
    }
    
    fn read_coordinates(&self) -> Result<CoordinateArrays, DataReaderError> {
        // Create placeholder coordinate arrays
        let lat_grid = Array2::<f32>::zeros((50, 50));
        let lon_grid = Array2::<f32>::zeros((50, 50));
        let levels = Array1::<f32>::from(vec![1000.0, 900.0, 800.0, 700.0, 600.0, 500.0, 400.0, 300.0, 200.0, 100.0]);
        
        Ok((lat_grid, lon_grid, levels))
    }
    
    fn get_global_attributes(&self) -> Result<Attributes, DataReaderError> {
        let mut attrs = Attributes::new();
        attrs.insert("title".to_string(), AttributeValue::String("NetCDF Test Data".to_string()));
        attrs.insert("source".to_string(), AttributeValue::String("QIBT Rust".to_string()));
        Ok(attrs)
    }
    
    fn get_variable_attributes(&self, variable_name: &str) -> Result<Attributes, DataReaderError> {
        let mut attrs = Attributes::new();
        
        match variable_name {
            "U" | "V" | "W" => {
                attrs.insert("units".to_string(), AttributeValue::String("m/s".to_string()));
                attrs.insert("_FillValue".to_string(), AttributeValue::Float(-9999.0));
            },
            "T" => {
                attrs.insert("units".to_string(), AttributeValue::String("K".to_string()));
                attrs.insert("_FillValue".to_string(), AttributeValue::Float(-9999.0));
            },
            "QVAPOR" => {
                attrs.insert("units".to_string(), AttributeValue::String("kg/kg".to_string()));
                attrs.insert("_FillValue".to_string(), AttributeValue::Float(-9999.0));
            },
            _ => {}
        }
        
        Ok(attrs)
    }
    
    fn get_global_attribute(&self, attribute_name: &str) -> Result<AttributeValue, DataReaderError> {
        match attribute_name {
            "title" => Ok(AttributeValue::String("NetCDF Test Data".to_string())),
            "source" => Ok(AttributeValue::String("QIBT Rust".to_string())),
            _ => Err(DataReaderError::MissingAttribute(attribute_name.to_string())),
        }
    }
    
    fn get_variable_attribute(
        &self,
        variable_name: &str,
        attribute_name: &str
    ) -> Result<AttributeValue, DataReaderError> {
        match (variable_name, attribute_name) {
            ("U" | "V" | "W", "units") => Ok(AttributeValue::String("m/s".to_string())),
            ("T", "units") => Ok(AttributeValue::String("K".to_string())),
            ("QVAPOR", "units") => Ok(AttributeValue::String("kg/kg".to_string())),
            (_, "_FillValue") => Ok(AttributeValue::Float(-9999.0)),
            _ => Err(DataReaderError::MissingAttribute(format!("{}.{}", variable_name, attribute_name))),
        }
    }
    
    fn get_path(&self) -> Option<PathBuf> {
        Some(PathBuf::from(&self.file_path))
    }
}
