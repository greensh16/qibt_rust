use crate::trajectory::TrajectoryPoint;
use std::path::Path;
use std::collections::HashMap;

/// Error types for data writing operations
#[derive(Debug, Clone)]
pub enum WriteError {
    IoError(String),
    FormatError(String),
    UnsupportedFeature(String),
    InvalidData(String),
}

impl std::fmt::Display for WriteError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WriteError::IoError(msg) => write!(f, "IO error: {}", msg),
            WriteError::FormatError(msg) => write!(f, "Format error: {}", msg),
            WriteError::UnsupportedFeature(msg) => write!(f, "Unsupported feature: {}", msg),
            WriteError::InvalidData(msg) => write!(f, "Invalid data: {}", msg),
        }
    }
}

impl std::error::Error for WriteError {}

/// Metadata for trajectory output files
#[derive(Debug, Clone)]
pub struct TrajectoryMetadata {
    pub start_longitude: f64,
    pub start_latitude: f64,
    pub start_pressure: f64,
    pub start_time: f64,
    pub integration_time_step: f64,
    pub total_integration_time: f64,
    pub meteorological_data_source: String,
    pub creation_time: String,
    pub global_attributes: HashMap<String, String>,
}

impl Default for TrajectoryMetadata {
    fn default() -> Self {
        Self {
            start_longitude: 0.0,
            start_latitude: 0.0,
            start_pressure: 101325.0,
            start_time: 0.0,
            integration_time_step: 600.0,
            total_integration_time: 86400.0,
            meteorological_data_source: "unknown".to_string(),
            creation_time: chrono::Utc::now().format("%Y-%m-%dT%H:%M:%S%.3fZ").to_string(),
            global_attributes: HashMap::new(),
        }
    }
}

/// Generic trait for writing trajectory data to different formats
pub trait DataWriter {
    /// Create or initialize the output file/dataset
    fn create(&mut self, expected_trajectories: usize, expected_time_steps: usize) -> Result<(), WriteError>;
    
    /// Write a single trajectory
    fn write_trajectory(&mut self, trajectory_id: u32, trajectory: &[TrajectoryPoint]) -> Result<(), WriteError>;
    
    /// Write multiple trajectories at once
    fn write_trajectories(&mut self, trajectories: &[(u32, Vec<TrajectoryPoint>)]) -> Result<(), WriteError>;
    
    /// Set global metadata and attributes
    fn set_metadata(&mut self, metadata: &TrajectoryMetadata) -> Result<(), WriteError>;
    
    /// Add custom global attribute
    fn add_global_attribute(&mut self, name: &str, value: &str) -> Result<(), WriteError>;
    
    /// Finalize and close the output file
    fn close(&mut self) -> Result<(), WriteError>;
    
    /// Get the output file path
    fn get_output_path(&self) -> &str;
}

/// Output format enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    NetCdf,
    Zarr,
    Ascii,
}

impl OutputFormat {
    /// Detect output format from file extension
    pub fn from_path(path: &Path) -> Self {
        match path.extension().and_then(|s| s.to_str()) {
            Some("nc") | Some("netcdf") | Some("nc4") => OutputFormat::NetCdf,
            Some("zarr") => OutputFormat::Zarr,
            Some("txt") | Some("ascii") | Some("csv") => OutputFormat::Ascii,
            _ => {
                // Check if path looks like a directory (for Zarr)
                if path.to_string_lossy().contains(".zarr") || 
                   (path.exists() && path.is_dir()) {
                    OutputFormat::Zarr
                } else {
                    OutputFormat::NetCdf // Default
                }
            }
        }
    }
    
    /// Get file extension for this format
    pub fn extension(&self) -> &'static str {
        match self {
            OutputFormat::NetCdf => "nc",
            OutputFormat::Zarr => "zarr",
            OutputFormat::Ascii => "txt",
        }
    }
}

impl std::fmt::Display for OutputFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OutputFormat::NetCdf => write!(f, "netcdf"),
            OutputFormat::Zarr => write!(f, "zarr"),
            OutputFormat::Ascii => write!(f, "ascii"),
        }
    }
}

impl std::str::FromStr for OutputFormat {
    type Err = String;
    
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "netcdf" | "nc" | "netcdf4" => Ok(OutputFormat::NetCdf),
            "zarr" => Ok(OutputFormat::Zarr),
            "ascii" | "txt" | "csv" => Ok(OutputFormat::Ascii),
            _ => Err(format!("Unknown output format: {}", s))
        }
    }
}

/// Factory function to create appropriate writer for the format
pub fn create_writer(output_path: &Path, format: OutputFormat) -> Result<Box<dyn DataWriter>, WriteError> {
    match format {
        OutputFormat::NetCdf => {
            let writer = crate::data_io::writer::NetCDFTrajectoryWriter::new(output_path)?;
            Ok(Box::new(writer))
        },
        OutputFormat::Zarr => {
            #[cfg(feature = "zarr")]
            {
                let writer = crate::data_io::zarr_writer::ZarrTrajectoryWriter::new(output_path)?;
                Ok(Box::new(writer))
            }
            #[cfg(not(feature = "zarr"))]
            {
                Err(WriteError::UnsupportedFeature("Zarr output requires 'zarr' feature".to_string()))
            }
        },
        OutputFormat::Ascii => {
            let writer = crate::data_io::ascii_writer::AsciiTrajectoryWriter::new(output_path)?;
            Ok(Box::new(writer))
        },
    }
}

/// Convenience function to auto-detect format and create writer
pub fn create_writer_auto(output_path: &Path) -> Result<Box<dyn DataWriter>, WriteError> {
    let format = OutputFormat::from_path(output_path);
    create_writer(output_path, format)
}
