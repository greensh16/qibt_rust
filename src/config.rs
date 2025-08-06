use chrono::{DateTime, NaiveDateTime, Utc};
use clap::{Arg, Command};
use std::path::PathBuf;
use std::sync::Arc;
use std::collections::HashMap;

/// Cloud authentication configuration
#[cfg(feature = "zarr")]
#[derive(Debug, Clone, Default)]
pub struct CloudAuth {
    /// AWS access key ID
    pub aws_access_key_id: Option<String>,
    /// AWS secret access key
    pub aws_secret_access_key: Option<String>,
    /// AWS session token (for temporary credentials)
    pub aws_session_token: Option<String>,
    /// AWS region
    pub aws_region: Option<String>,
    /// GCP service account key path
    pub gcp_service_account_key: Option<String>,
    /// GCP project ID
    pub gcp_project_id: Option<String>,
    /// Azure storage account name
    pub azure_storage_account: Option<String>,
    /// Azure storage account key
    pub azure_storage_key: Option<String>,
    /// Custom endpoint URL (for S3-compatible services)
    pub endpoint_url: Option<String>,
    /// Additional authentication headers
    pub custom_headers: HashMap<String, String>,
}

#[cfg(feature = "zarr")]
#[derive(Debug, Clone)]
pub struct StreamingConfig {
    /// Chunk size for streaming reads (in bytes)
    pub chunk_size: usize,
    /// Maximum number of concurrent requests
    pub max_concurrent_requests: usize,
    /// Request timeout in seconds
    pub request_timeout_secs: u64,
    /// Maximum number of retries for failed requests
    pub max_retries: usize,
    /// Base delay between retries in milliseconds
    pub retry_delay_ms: u64,
    /// Whether to use exponential backoff for retries
    pub exponential_backoff: bool,
    /// Buffer size for in-memory caching
    pub buffer_size: usize,
    /// Whether to enable range request optimization
    pub enable_range_requests: bool,
}

#[cfg(feature = "zarr")]
impl Default for StreamingConfig {
    fn default() -> Self {
        Self {
            chunk_size: 4 * 1024 * 1024, // 4MB chunks
            max_concurrent_requests: 10,
            request_timeout_secs: 30,
            max_retries: 3,
            retry_delay_ms: 100,
            exponential_backoff: true,
            buffer_size: 64 * 1024 * 1024, // 64MB buffer
            enable_range_requests: true,
        }
    }
}

/// Global configuration constants - porting MODULE global_data values
#[derive(Clone, Debug)]
pub struct Constants {
    /// Gravitational acceleration (m/s²) - from Fortran: g = 9.8
    pub g: f64,
    /// Gas constant for dry air (J/(kg·K)) - from Fortran: Rd = 287.053
    pub r_dry: f64,
    /// Specific heat at constant pressure (J/(kg·K)) - from Fortran: Cp = 1004.67
    pub cp: f64,
    /// Reference pressure (Pa) - from Fortran: P0 = 100000
    pub p0: f64,
    /// Earth's radius (m)
    pub earth_radius: f64,

    // Thermodynamic constants from Fortran
    /// Latent heat of vaporization of water (J/kg) - from Fortran: Lv = 2.25E6
    pub lv: f64,
    /// Gas constant of water vapor (J/(kg·K)) - from Fortran: Rv = 461.5
    pub rv: f64,
    /// Heat capacity of liquid water at ~-20C (J/(kg·K)) - from Fortran: Cl = 4400
    pub cl: f64,
    /// Average distance of 1 degree lat (km) - from Fortran: deg_dist = 111
    pub deg_dist: f64,
    /// Pi constant - from Fortran: pi = 3.14159265
    pub pi: f64,

    // Simulation parameters from MODULE global_data
    /// Time step in minutes (TSTEP equivalent) - from Fortran: tstep = 10
    pub tstep_minutes: i32,
    /// Number of parcels for ensemble runs (N_PARCELS equivalent) - from Fortran: nparcels = 100
    pub n_parcels: usize,
    /// Number of back-tracking days - from Fortran: totbtadays = 30
    pub tot_bta_days: i32,
    /// Minimum daily precipitation to process (mm) - from Fortran: minpre = 2
    pub min_precip: f64,
    /// Boundary layers to ignore - from Fortran: bdy = 6
    pub boundary_layers: i32,
    /// Number of parallel threads - from Fortran: numthreads = 8
    pub num_threads: usize,
    /// Minimum change in parcel mixing ratio (kg/kg) - from Fortran: min_del_q = 0.0001
    pub min_del_q: f64,

    // Data input parameters
    /// Input data time step in minutes - from Fortran: datatstep = 180
    pub data_timestep_minutes: i32,

    // Model domain boundaries
    pub lon_min: f64,
    pub lon_max: f64,
    pub lat_min: f64,
    pub lat_max: f64,
    /// Vertical domain boundaries (Pa)
    pub pressure_top: f64,
    pub pressure_surface: f64,
    /// Maximum integration time (hours)
    pub max_integration_time: f64,

    // Quality control thresholds
    pub max_horizontal_speed: f64, // m/s
    pub max_vertical_speed: f64,   // m/s
    pub min_pressure: f64,         // Pa
    pub max_pressure: f64,         // Pa

    // Flags and options
    /// Process storm peaks vs whole days - from Fortran: peak = .FALSE.
    pub process_peaks: bool,
    /// Only calculate trajectories for watershed - from Fortran: wshed = .TRUE.
    pub watershed_only: bool,
}

impl Default for Constants {
    fn default() -> Self {
        Self {
            // Physical constants
            g: 9.81,
            r_dry: 287.0,
            cp: 1004.0,
            p0: 100000.0,
            earth_radius: 6371000.0,

            // Thermodynamic constants from Fortran
            lv: 2.25e6,               // Latent heat of vaporization (J/kg)
            rv: 461.5,                // Gas constant of water vapor (J/(kg·K))
            cl: 4400.0,               // Heat capacity of liquid water (J/(kg·K))
            deg_dist: 111.0,          // Average distance of 1 degree lat (km)
            pi: std::f64::consts::PI, // Pi constant

            // Simulation parameters (porting from MODULE global_data)
            tstep_minutes: 10,  // Time step in minutes
            n_parcels: 100,     // Default number of parcels
            tot_bta_days: 30,   // Number of back-tracking days
            min_precip: 2.0,    // Minimum daily precipitation (mm)
            boundary_layers: 6, // Boundary layers to ignore
            num_threads: 8,     // Number of parallel threads
            min_del_q: 0.0001,  // Minimum change in parcel mixing ratio (kg/kg)

            // Data input parameters
            data_timestep_minutes: 180, // Input data time step in minutes

            // Model domain (global by default)
            lon_min: -180.0,
            lon_max: 180.0,
            lat_min: -90.0,
            lat_max: 90.0,

            // Vertical domain (surface to ~10 hPa)
            pressure_top: 1000.0,       // 10 hPa
            pressure_surface: 101325.0, // 1013.25 hPa

            // Integration limits
            max_integration_time: 240.0, // 10 days

            // Quality control
            max_horizontal_speed: 200.0, // 200 m/s (~720 km/h)
            max_vertical_speed: 10.0,    // 10 m/s
            min_pressure: 100.0,         // 1 hPa minimum
            max_pressure: 110000.0,      // 1100 hPa maximum

            // Flags and options
            process_peaks: false, // Process storm peaks vs whole days
            watershed_only: true, // Only calculate trajectories for watershed
        }
    }
}

/// Main configuration structure with CLI support
#[derive(Clone, Debug)]
pub struct Config {
    /// Physical constants
    pub constants: Constants,

    // Input/Output paths (mutable via CLI)
    /// Input mode (multifile or singlefile)
    pub input_mode: InputMode,
    /// Input meteorological data path (file or directory)
    pub input_path: PathBuf,
    /// Output directory for trajectory files
    pub output_dir: PathBuf,
    /// Output trajectory file name pattern
    pub output_pattern: String,

    // Time configuration (mutable via CLI)
    /// Starting date and time
    pub start_date: DateTime<Utc>,
    /// Ending date and time
    pub end_date: DateTime<Utc>,
    /// Starting time (Julian day)
    pub start_time: f64,
    /// Ending time (Julian day)
    pub end_time: f64,

    // Release configuration
    /// Starting latitude (degrees)
    pub start_lat: f64,
    /// Starting longitude (degrees)
    pub start_lon: f64,
    /// Starting pressure level (Pa)
    pub start_pressure: f64,

    // Simulation parameters
    /// Trajectory length (hours)
    pub trajectory_length: f64,
    /// Time step (seconds) - defaults to constants.tstep
    pub time_step: f64,
    /// Number of parcels for ensemble runs
    pub num_parcels: usize,
    /// Number of parallel threads
    pub num_threads: usize,
    /// Cloud authentication configuration
    #[cfg(feature = "zarr")]
    pub cloud_auth: CloudAuth,
    /// Streaming configuration for remote data access
    #[cfg(feature = "zarr")]
    pub streaming_config: StreamingConfig,

    // Model configuration
    /// Integration method
    pub integration_method: IntegrationMethod,
    /// Output frequency (hours)
    pub output_frequency: f64,
    /// Enable quality control checks
    pub enable_qc: bool,
    /// Verbose output
    pub verbose: bool,
}

/// Available input modes
#[derive(Clone, Debug)]
pub enum InputMode {
    /// Process multiple files from a directory
    MultiFile,
    /// Process a single file
    SingleFile,
}

/// Available integration methods
#[derive(Clone, Debug)]
pub enum IntegrationMethod {
    /// Simple backward Euler
    BackwardEuler,
    /// Implicit backward integration
    ImplicitBackward,
    /// Fourth-order Runge-Kutta
    RungeKutta4,
}

impl Default for Config {
    fn default() -> Self {
        let constants = Constants::default();
        let default_start = Utc::now();
        let default_end = default_start + chrono::Duration::hours(120); // 5 days

        Self {
            constants: constants.clone(),
            input_mode: InputMode::MultiFile,
            input_path: PathBuf::from("./input"),
            output_dir: PathBuf::from("./output"),
            output_pattern: String::from("trajectory_%Y%m%d_%H%M.nc"),

            start_date: default_start,
            end_date: default_end,
            start_time: Self::datetime_to_julian(&default_start),
            end_time: Self::datetime_to_julian(&default_end),

            start_lat: 45.0,
            start_lon: -100.0,
            start_pressure: 50000.0, // 500 hPa

            trajectory_length: 120.0, // 5 days
            time_step: 3600.0,        // 1 hour default
            num_parcels: constants.n_parcels,
            num_threads: 4,
            #[cfg(feature = "zarr")]
            cloud_auth: CloudAuth::default(),
            #[cfg(feature = "zarr")]
            streaming_config: StreamingConfig::default(),

            integration_method: IntegrationMethod::ImplicitBackward,
            output_frequency: 1.0, // hourly output
            enable_qc: true,
            verbose: false,
        }
    }
}

impl Config {
    /// Parse configuration from command line arguments
    pub fn from_args() -> Result<Self, String> {
        let app = Command::new("qibt_rust")
            .version("0.1.0")
            .author("Back-trajectory analysis tool")
            .about("Quasi-isentropic back-trajectory model")
            .arg(
                Arg::new("input-path")
                    .short('i')
                    .long("input-path")
                    .value_name("PATH")
                    .help("Input meteorological data path (file or directory)")
                    .required(true),
            )
            .arg(
                Arg::new("input-mode")
                    .long("input-mode")
                    .value_name("MODE")
                    .help("Input mode: multifile (directory) or singlefile")
                    .value_parser(["multifile", "singlefile"])
                    .default_value("multifile"),
            )
            .arg(
                Arg::new("output-dir")
                    .short('o')
                    .long("output-dir")
                    .value_name("DIR")
                    .help("Output directory for trajectory files")
                    .default_value("./output"),
            )
            .arg(
                Arg::new("output-pattern")
                    .short('p')
                    .long("output-pattern")
                    .value_name("PATTERN")
                    .help("Output file name pattern (strftime format)")
                    .default_value("trajectory_%Y%m%d_%H%M.nc"),
            )
            .arg(
                Arg::new("start-date")
                    .short('s')
                    .long("start-date")
                    .value_name("DATETIME")
                    .help("Start date and time (YYYY-MM-DD HH:MM:SS)")
                    .required(true),
            )
            .arg(
                Arg::new("end-date")
                    .short('e')
                    .long("end-date")
                    .value_name("DATETIME")
                    .help("End date and time (YYYY-MM-DD HH:MM:SS)")
                    .required(true),
            )
            .arg(
                Arg::new("start-lat")
                    .long("start-lat")
                    .value_name("DEGREES")
                    .help("Starting latitude (degrees)")
                    .default_value("45.0"),
            )
            .arg(
                Arg::new("start-lon")
                    .long("start-lon")
                    .value_name("DEGREES")
                    .help("Starting longitude (degrees)")
                    .default_value("-100.0"),
            )
            .arg(
                Arg::new("start-pressure")
                    .long("start-pressure")
                    .value_name("PA")
                    .help("Starting pressure level (Pa)")
                    .default_value("50000"),
            )
            .arg(
                Arg::new("trajectory-length")
                    .short('l')
                    .long("trajectory-length")
                    .value_name("HOURS")
                    .help("Trajectory length (hours)")
                    .default_value("120"),
            )
            .arg(
                Arg::new("time-step")
                    .short('t')
                    .long("time-step")
                    .value_name("SECONDS")
                    .help("Time step (seconds)")
                    .default_value("3600"),
            )
            .arg(
                Arg::new("num-parcels")
                    .short('n')
                    .long("num-parcels")
                    .value_name("COUNT")
                    .help("Number of parcels for ensemble runs")
                    .default_value("1000"),
            )
            .arg(
                Arg::new("num-threads")
                    .short('j')
                    .long("num-threads")
                    .value_name("COUNT")
                    .help("Number of parallel threads")
                    .default_value("4"),
            )
            .arg(
                Arg::new("integration-method")
                    .short('m')
                    .long("integration-method")
                    .value_name("METHOD")
                    .help("Integration method")
                    .value_parser(["euler", "implicit", "rk4"])
                    .default_value("implicit"),
            )
            .arg(
                Arg::new("output-frequency")
                    .short('f')
                    .long("output-frequency")
                    .value_name("HOURS")
                    .help("Output frequency (hours)")
                    .default_value("1.0"),
            )
            .arg(
                Arg::new("no-qc")
                    .long("no-qc")
                    .help("Disable quality control checks")
                    .action(clap::ArgAction::SetTrue),
            )
            .arg(
                Arg::new("verbose")
                    .short('v')
                    .long("verbose")
                    .help("Enable verbose output")
                    .action(clap::ArgAction::SetTrue),
            );

        let matches = app.try_get_matches().map_err(|e| e.to_string())?;

        // Parse required arguments
        let input_path = PathBuf::from(matches.get_one::<String>("input-path").unwrap());
        let output_dir = PathBuf::from(matches.get_one::<String>("output-dir").unwrap());
        let output_pattern = matches.get_one::<String>("output-pattern").unwrap().clone();
        
        // Parse input mode
        let input_mode = match matches.get_one::<String>("input-mode").unwrap().as_str() {
            "multifile" => InputMode::MultiFile,
            "singlefile" => InputMode::SingleFile,
            _ => return Err("Invalid input mode".to_string()),
        };

        // Parse dates
        let start_date = Self::parse_datetime(matches.get_one::<String>("start-date").unwrap())?;
        let end_date = Self::parse_datetime(matches.get_one::<String>("end-date").unwrap())?;

        // Parse numeric parameters
        let start_lat: f64 = matches
            .get_one::<String>("start-lat")
            .unwrap()
            .parse()
            .map_err(|_| "Invalid start latitude")?;
        let start_lon: f64 = matches
            .get_one::<String>("start-lon")
            .unwrap()
            .parse()
            .map_err(|_| "Invalid start longitude")?;
        let start_pressure: f64 = matches
            .get_one::<String>("start-pressure")
            .unwrap()
            .parse()
            .map_err(|_| "Invalid start pressure")?;
        let trajectory_length: f64 = matches
            .get_one::<String>("trajectory-length")
            .unwrap()
            .parse()
            .map_err(|_| "Invalid trajectory length")?;
        let time_step: f64 = matches
            .get_one::<String>("time-step")
            .unwrap()
            .parse()
            .map_err(|_| "Invalid time step")?;
        let num_parcels: usize = matches
            .get_one::<String>("num-parcels")
            .unwrap()
            .parse()
            .map_err(|_| "Invalid number of parcels")?;
        let num_threads: usize = matches
            .get_one::<String>("num-threads")
            .unwrap()
            .parse()
            .map_err(|_| "Invalid number of threads")?;
        let output_frequency: f64 = matches
            .get_one::<String>("output-frequency")
            .unwrap()
            .parse()
            .map_err(|_| "Invalid output frequency")?;

        // Parse integration method
        let integration_method = match matches
            .get_one::<String>("integration-method")
            .unwrap()
            .as_str()
        {
            "euler" => IntegrationMethod::BackwardEuler,
            "implicit" => IntegrationMethod::ImplicitBackward,
            "rk4" => IntegrationMethod::RungeKutta4,
            _ => return Err("Invalid integration method".to_string()),
        };

        // Parse flags
        let enable_qc = !matches.get_flag("no-qc");
        let verbose = matches.get_flag("verbose");

        let config = Self {
            constants: Constants::default(),
            input_mode,
            input_path,
            output_dir,
            output_pattern,
            start_date,
            end_date,
            start_time: Self::datetime_to_julian(&start_date),
            end_time: Self::datetime_to_julian(&end_date),
            start_lat,
            start_lon,
            start_pressure,
            trajectory_length,
            time_step,
            num_parcels,
            num_threads,
            #[cfg(feature = "zarr")]
            cloud_auth: CloudAuth::default(),
            #[cfg(feature = "zarr")]
            streaming_config: StreamingConfig::default(),
            integration_method,
            output_frequency,
            enable_qc,
            verbose,
        };

        // Validate configuration
        config.validate()?;

        Ok(config)
    }

    /// Create an Arc<Config> for thread-safe sharing
    pub fn into_arc(self) -> Arc<Self> {
        Arc::new(self)
    }

    /// Parse datetime string in format "YYYY-MM-DD HH:MM:SS"
    fn parse_datetime(datetime_str: &str) -> Result<DateTime<Utc>, String> {
        NaiveDateTime::parse_from_str(datetime_str, "%Y-%m-%d %H:%M:%S")
            .map_err(|_| {
                format!(
                    "Invalid datetime format: {}. Expected: YYYY-MM-DD HH:MM:SS",
                    datetime_str
                )
            })
            .map(|dt| DateTime::<Utc>::from_naive_utc_and_offset(dt, Utc))
    }

    /// Convert DateTime to Julian Day Number
    fn datetime_to_julian(dt: &DateTime<Utc>) -> f64 {
        // Julian Day Number calculation
        let timestamp = dt.timestamp() as f64;
        // Unix epoch (1970-01-01 00:00:00 UTC) corresponds to Julian Day 2440587.5
        2440587.5 + timestamp / 86400.0
    }

    /// Convert Julian Day Number to DateTime
    pub fn julian_to_datetime(jd: f64) -> DateTime<Utc> {
        let timestamp = (jd - 2440587.5) * 86400.0;
        DateTime::<Utc>::from_timestamp(timestamp as i64, 0).unwrap_or(Utc::now())
    }

    /// Create a Config for testing purposes (bypasses CLI parsing)
    #[cfg(test)]
    pub fn for_testing(
        input_mode: InputMode,
        input_path: std::path::PathBuf,
    ) -> Result<Self, String> {
        let config = Self {
            input_mode,
            input_path,
            ..Self::default()
        };
        config.validate()?;
        Ok(config)
    }
    
    /// Validate configuration parameters
    pub fn validate(&self) -> Result<(), String> {
        // Validate coordinate ranges
        if self.start_lat < -90.0 || self.start_lat > 90.0 {
            return Err("Latitude must be between -90 and 90 degrees".to_string());
        }
        if self.start_lon < -180.0 || self.start_lon > 180.0 {
            return Err("Longitude must be between -180 and 180 degrees".to_string());
        }
        if self.start_pressure <= 0.0 {
            return Err("Pressure must be positive".to_string());
        }
        if self.time_step <= 0.0 {
            return Err("Time step must be positive".to_string());
        }
        
        // Validate input path exists
        if !self.input_path.exists() {
            return Err(format!("Input path does not exist: {}", self.input_path.display()));
        }
        
        // Validate input path is compatible with input mode
        match self.input_mode {
            InputMode::MultiFile => {
                if !self.input_path.is_dir() {
                    return Err(format!(
                        "Input mode is multifile but path is not a directory: {}", 
                        self.input_path.display()
                    ));
                }
            }
            InputMode::SingleFile => {
                if !self.input_path.is_file() {
                    return Err(format!(
                        "Input mode is singlefile but path is not a file: {}", 
                        self.input_path.display()
                    ));
                }
            }
        }
        
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;

    #[test]
    fn test_input_mode_default() {
        let config = Config::default();
        matches!(config.input_mode, InputMode::MultiFile);
    }

    #[test]
    fn test_config_multifile_validation_success() {
        // Create a temporary directory for testing
        let test_dir = "test_temp_dir";
        fs::create_dir_all(test_dir).unwrap();
        
        let result = Config::for_testing(
            InputMode::MultiFile,
            Path::new(test_dir).to_path_buf(),
        );
        
        assert!(result.is_ok());
        let config = result.unwrap();
        matches!(config.input_mode, InputMode::MultiFile);
        
        // Clean up
        fs::remove_dir_all(test_dir).unwrap();
    }

    #[test]
    fn test_config_singlefile_validation_success() {
        // Use Cargo.toml which should exist
        let result = Config::for_testing(
            InputMode::SingleFile,
            Path::new("Cargo.toml").to_path_buf(),
        );
        
        assert!(result.is_ok());
        let config = result.unwrap();
        matches!(config.input_mode, InputMode::SingleFile);
    }

    #[test]
    fn test_config_multifile_validation_failure_with_file() {
        // Try to use a file with multifile mode
        let result = Config::for_testing(
            InputMode::MultiFile,
            Path::new("Cargo.toml").to_path_buf(),
        );
        
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Input mode is multifile but path is not a directory"));
    }

    #[test]
    fn test_config_singlefile_validation_failure_with_directory() {
        // Try to use a directory with singlefile mode
        let result = Config::for_testing(
            InputMode::SingleFile,
            Path::new("src").to_path_buf(),
        );
        
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Input mode is singlefile but path is not a file"));
    }

    #[test]
    fn test_config_validation_nonexistent_path() {
        let result = Config::for_testing(
            InputMode::MultiFile,
            Path::new("nonexistent_path_12345").to_path_buf(),
        );
        
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Input path does not exist"));
    }

    #[test]
    fn test_input_mode_debug() {
        let multifile = InputMode::MultiFile;
        let singlefile = InputMode::SingleFile;
        
        assert_eq!(format!("{:?}", multifile), "MultiFile");
        assert_eq!(format!("{:?}", singlefile), "SingleFile");
    }
}
