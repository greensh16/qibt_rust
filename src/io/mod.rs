use ndarray::{Array1, Array2, Array4};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{Read, BufReader};
use thiserror::Error;
use chrono::{DateTime, Utc};

#[cfg(feature = "zarr")]
mod zarr_reader;
// #[cfg(feature = "zarr")]
// mod cloud_auth;
// #[cfg(feature = "zarr")]
// mod streaming_reader;

/// Generic error type for data readers
#[derive(Error, Debug)]
pub enum DataReaderError {
    #[error("NetCDF error: {0}")]
    Netcdf(String),
    
    #[error("Zarr error: {0}")]
    Zarr(String),
    
    #[error("Variable not found: {0}")]
    MissingVariable(String),
    
    #[error("Dimension not found: {0}")]
    MissingDimension(String),
    
    #[error("Attribute not found: {0}")]
    MissingAttribute(String),
    
    #[error("Data conversion error: {0}")]
    ConversionError(String),
    
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    
    #[error("Invalid file format: {0}")]
    InvalidFormat(String),
    
    #[error("File not found: {0}")]
    FileNotFound(String),
    
    #[error("Unsupported operation: {0}")]
    UnsupportedOperation(String),
}

/// Metadata about a variable in the dataset
#[derive(Debug, Clone)]
pub struct VariableInfo {
    /// Variable name
    pub name: String,
    /// Variable dimensions
    pub dimensions: Vec<String>,
    /// Variable shape
    pub shape: Vec<usize>,
    /// Data type
    pub dtype: String,
    /// Units attribute if available
    pub units: Option<String>,
    /// Long name or description if available
    pub long_name: Option<String>,
}

/// Metadata about a dimension in the dataset
#[derive(Debug, Clone)]
pub struct DimensionInfo {
    /// Dimension name
    pub name: String,
    /// Dimension size
    pub size: usize,
    /// Whether this is an unlimited dimension
    pub is_unlimited: bool,
}

/// Global and variable attributes
pub type Attributes = HashMap<String, AttributeValue>;

/// Type alias for coordinate arrays (lat, lon, levels)
pub type CoordinateArrays = (Array2<f32>, Array2<f32>, Array1<f32>);

/// Supported attribute value types
#[derive(Debug, Clone)]
pub enum AttributeValue {
    String(String),
    Int(i32),
    Float(f32),
    Double(f64),
    IntArray(Vec<i32>),
    FloatArray(Vec<f32>),
    DoubleArray(Vec<f64>),
}

/// File metadata structure
#[derive(Debug, Clone)]
pub struct FileMetadata {
    /// Available dimensions
    pub dimensions: Vec<DimensionInfo>,
    /// Available variables
    pub variables: Vec<VariableInfo>,
    /// Global attributes
    pub global_attributes: Attributes,
}

/// Generic trait for reading scientific data files (NetCDF, Zarr, etc.)
/// 
/// This trait provides a common interface for reading multidimensional array data
/// from different file formats, with methods for listing variables, reading data slices,
/// and accessing metadata.
pub trait DataReader {
    /// Open a file or dataset for reading
    /// 
    /// # Arguments
    /// * `path` - Path to the file or dataset
    /// 
    /// # Returns
    /// * `Result<(), DataReaderError>` - Success or error
    fn open(&mut self, path: &Path) -> Result<(), DataReaderError>;
    
    /// Check if the reader is currently open and ready to read data
    fn is_open(&self) -> bool;
    
    /// Close the current file/dataset
    fn close(&mut self);
    
    /// List all available variables in the dataset
    /// 
    /// # Returns
    /// * `Result<Vec<String>, DataReaderError>` - Variable names or error
    fn list_variables(&self) -> Result<Vec<String>, DataReaderError>;
    
    /// Get detailed information about a specific variable
    /// 
    /// # Arguments
    /// * `variable_name` - Name of the variable
    /// 
    /// # Returns
    /// * `Result<VariableInfo, DataReaderError>` - Variable metadata or error
    fn get_variable_info(&self, variable_name: &str) -> Result<VariableInfo, DataReaderError>;
    
    /// List all dimensions in the dataset
    /// 
    /// # Returns
    /// * `Result<Vec<String>, DataReaderError>` - Dimension names or error  
    fn list_dimensions(&self) -> Result<Vec<String>, DataReaderError>;
    
    /// Get information about a specific dimension
    /// 
    /// # Arguments
    /// * `dimension_name` - Name of the dimension
    /// 
    /// # Returns
    /// * `Result<DimensionInfo, DataReaderError>` - Dimension metadata or error
    fn get_dimension_info(&self, dimension_name: &str) -> Result<DimensionInfo, DataReaderError>;
    
    /// Get complete file metadata including all variables and dimensions
    /// 
    /// # Returns
    /// * `Result<FileMetadata, DataReaderError>` - Complete metadata or error
    fn get_metadata(&self) -> Result<FileMetadata, DataReaderError>;
    
    /// Read a complete variable as a 4D array
    /// 
    /// This method reads the entire variable and returns it as a 4D ndarray.
    /// For variables with fewer than 4 dimensions, the array is padded with singleton dimensions.
    /// 
    /// # Arguments
    /// * `variable_name` - Name of the variable to read
    /// 
    /// # Returns
    /// * `Result<Array4<f32>, DataReaderError>` - 4D array or error
    fn read_variable(&self, variable_name: &str) -> Result<Array4<f32>, DataReaderError>;
    
    /// Read a slice of a variable as a 4D array
    /// 
    /// This method allows reading only a portion of a variable, which is useful
    /// for large datasets where only a subset of the data is needed.
    /// 
    /// # Arguments
    /// * `variable_name` - Name of the variable to read
    /// * `indices` - Slice indices for each dimension (start, end, step)
    /// 
    /// # Returns
    /// * `Result<Array4<f32>, DataReaderError>` - 4D array slice or error
    fn read_variable_slice(
        &self, 
        variable_name: &str, 
        indices: &[(usize, usize, usize)]
    ) -> Result<Array4<f32>, DataReaderError>;
    
    /// Read multiple variables efficiently
    /// 
    /// This method can be more efficient than reading variables individually
    /// as it can optimize file access patterns.
    /// 
    /// # Arguments
    /// * `variable_names` - Names of variables to read
    /// 
    /// # Returns
    /// * `Result<HashMap<String, Array4<f32>>, DataReaderError>` - Map of variable data or error
    fn read_variables(&self, variable_names: &[&str]) -> Result<HashMap<String, Array4<f32>>, DataReaderError>;
    
    /// Read coordinate variables (typically longitude, latitude, levels)
    /// 
    /// This is a convenience method for reading spatial coordinate information.
    /// Returns 2D arrays for lat/lon and 1D array for vertical levels.
    /// 
    /// # Returns
    /// * `Result<(Array2<f32>, Array2<f32>, Array1<f32>), DataReaderError>` - (lat, lon, levels) or error
    fn read_coordinates(&self) -> Result<CoordinateArrays, DataReaderError>;
    
    /// Get global attributes from the file
    /// 
    /// # Returns
    /// * `Result<Attributes, DataReaderError>` - Global attributes or error
    fn get_global_attributes(&self) -> Result<Attributes, DataReaderError>;
    
    /// Get attributes for a specific variable
    /// 
    /// # Arguments
    /// * `variable_name` - Name of the variable
    /// 
    /// # Returns
    /// * `Result<Attributes, DataReaderError>` - Variable attributes or error
    fn get_variable_attributes(&self, variable_name: &str) -> Result<Attributes, DataReaderError>;
    
    /// Get a specific global attribute value
    /// 
    /// # Arguments
    /// * `attribute_name` - Name of the attribute
    /// 
    /// # Returns
    /// * `Result<AttributeValue, DataReaderError>` - Attribute value or error
    fn get_global_attribute(&self, attribute_name: &str) -> Result<AttributeValue, DataReaderError>;
    
    /// Get a specific variable attribute value
    /// 
    /// # Arguments
    /// * `variable_name` - Name of the variable
    /// * `attribute_name` - Name of the attribute
    /// 
    /// # Returns
    /// * `Result<AttributeValue, DataReaderError>` - Attribute value or error
    fn get_variable_attribute(
        &self, 
        variable_name: &str, 
        attribute_name: &str
    ) -> Result<AttributeValue, DataReaderError>;
    
    /// Check if a variable exists in the dataset
    /// 
    /// # Arguments
    /// * `variable_name` - Name of the variable to check
    /// 
    /// # Returns
    /// * `bool` - True if variable exists
    fn has_variable(&self, variable_name: &str) -> bool {
        self.list_variables()
            .map(|vars| vars.contains(&variable_name.to_string()))
            .unwrap_or(false)
    }
    
    /// Check if a dimension exists in the dataset
    /// 
    /// # Arguments
    /// * `dimension_name` - Name of the dimension to check
    /// 
    /// # Returns
    /// * `bool` - True if dimension exists
    fn has_dimension(&self, dimension_name: &str) -> bool {
        self.list_dimensions()
            .map(|dims| dims.contains(&dimension_name.to_string()))
            .unwrap_or(false)
    }
    
    /// Get the file path or dataset identifier
    /// 
    /// # Returns
    /// * `Option<PathBuf>` - Path if available
    fn get_path(&self) -> Option<PathBuf>;
}

/// Trait for readers that support time-based data access
/// 
/// This trait extends DataReader for datasets that have a time dimension,
/// providing methods to read data for specific time steps.
pub trait TimeSeriesDataReader: DataReader {
    /// Read data for a specific datetime
    /// 
    /// # Arguments
    /// * `datetime` - The datetime to read data for
    /// * `variables` - Variables to read
    /// 
    /// # Returns
    /// * `Result<HashMap<String, Array4<f32>>, DataReaderError>` - Variable data or error
    fn read_at_time(
        &self,
        datetime: &DateTime<Utc>,
        variables: &[&str],
    ) -> Result<HashMap<String, Array4<f32>>, DataReaderError>;
    
    /// Get available time steps in the dataset
    /// 
    /// # Returns
    /// * `Result<Vec<DateTime<Utc>>, DataReaderError>` - Available times or error
    fn get_time_steps(&self) -> Result<Vec<DateTime<Utc>>, DataReaderError>;
    
    /// Find the closest time step to a given datetime
    /// 
    /// # Arguments
    /// * `target_time` - Target datetime
    /// 
    /// # Returns
    /// * `Result<DateTime<Utc>, DataReaderError>` - Closest available time or error
    fn find_closest_time(&self, target_time: &DateTime<Utc>) -> Result<DateTime<Utc>, DataReaderError>;
}

/// Convenience functions for creating readers
/// Create a NetCDF reader
/// 
/// # Arguments
/// * `path` - Path to NetCDF file
/// 
/// # Returns
/// * `Result<Box<dyn DataReader>, DataReaderError>` - NetCDF reader or error
pub fn create_netcdf_reader<P: AsRef<Path>>(_path: P) -> Result<Box<dyn DataReader>, DataReaderError> {
    // This would be implemented when NetCDF reader is available
    Err(DataReaderError::UnsupportedOperation(
        "NetCDF reader not yet implemented".to_string()
    ))
}

/// Create a Zarr reader
/// 
/// # Arguments
/// * `path` - Path to Zarr dataset
/// 
/// # Returns
/// * `Result<Box<dyn DataReader>, DataReaderError>` - Zarr reader or error
#[cfg(feature = "zarr")]
pub fn create_zarr_reader<P: AsRef<Path>>(path: P) -> Result<Box<dyn DataReader>, DataReaderError> {
    let mut reader = zarr_reader::ZarrReader::new();
    reader.open(path.as_ref())?;
    Ok(Box::new(reader))
}

#[cfg(not(feature = "zarr"))]
pub fn create_zarr_reader<P: AsRef<Path>>(_path: P) -> Result<Box<dyn DataReader>, DataReaderError> {
    Err(DataReaderError::UnsupportedOperation(
        "Zarr support not compiled in. Enable the 'zarr' feature.".to_string()
    ))
}

/// Dataset structure that provides a unified interface for opening different data formats
/// 
/// This struct wraps the auto-detection functionality and provides a convenient
/// interface similar to other scientific data libraries.
pub struct Dataset {
    reader: Box<dyn DataReader>,
}

impl Dataset {
    /// Open a dataset with automatic format detection
    /// 
    /// This method implements the format detection strategy as specified:
    /// 1. First tries NetCDF magic bytes detection
    /// 2. If that fails, checks for Zarr format indicators:
    ///    - `.zarr` suffix in path
    ///    - Directory with `.zgroup` or `.zarray` files
    ///    - JSON header with Zarr metadata
    /// 
    /// Existing NetCDF behavior remains unchanged - any file detected as NetCDF
    /// will be handled by the NetCDF reader.
    /// 
    /// # Arguments
    /// * `path` - Path to data file/dataset (file path or URL)
    /// 
    /// # Returns
    /// * `Result<Dataset, DataReaderError>` - Dataset instance with appropriate reader
    /// 
    /// # Examples
    /// 
    /// ```rust
    /// use qibt_rust::io::Dataset;
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// // Opens NetCDF file (detected by magic bytes)
    /// // let netcdf_dataset = Dataset::open("/path/to/data.nc")?;
    /// 
    /// // Opens Zarr dataset (detected by directory structure)
    /// // let zarr_dataset = Dataset::open("/path/to/data.zarr")?;
    /// 
    /// // Opens Zarr dataset (detected by .zgroup file)
    /// // let zarr_dir = Dataset::open("/path/to/zarr_directory")?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, DataReaderError> {
        let reader = create_reader(path)?;
        Ok(Dataset { reader })
    }
    
    /// Get the underlying DataReader
    /// 
    /// This provides access to the full DataReader interface for advanced usage.
    /// 
    /// # Returns
    /// * `&dyn DataReader` - Reference to the underlying reader
    pub fn reader(&self) -> &dyn DataReader {
        self.reader.as_ref()
    }
    
    /// Get a mutable reference to the underlying DataReader
    /// 
    /// # Returns
    /// * `&mut dyn DataReader` - Mutable reference to the underlying reader
    pub fn reader_mut(&mut self) -> &mut dyn DataReader {
        self.reader.as_mut()
    }
}

// Implement DataReader trait for Dataset to provide direct access to reader methods
impl DataReader for Dataset {
    fn open(&mut self, path: &Path) -> Result<(), DataReaderError> {
        self.reader.open(path)
    }
    
    fn is_open(&self) -> bool {
        self.reader.is_open()
    }
    
    fn close(&mut self) {
        self.reader.close()
    }
    
    fn list_variables(&self) -> Result<Vec<String>, DataReaderError> {
        self.reader.list_variables()
    }
    
    fn get_variable_info(&self, variable_name: &str) -> Result<VariableInfo, DataReaderError> {
        self.reader.get_variable_info(variable_name)
    }
    
    fn list_dimensions(&self) -> Result<Vec<String>, DataReaderError> {
        self.reader.list_dimensions()
    }
    
    fn get_dimension_info(&self, dimension_name: &str) -> Result<DimensionInfo, DataReaderError> {
        self.reader.get_dimension_info(dimension_name)
    }
    
    fn get_metadata(&self) -> Result<FileMetadata, DataReaderError> {
        self.reader.get_metadata()
    }
    
    fn read_variable(&self, variable_name: &str) -> Result<Array4<f32>, DataReaderError> {
        self.reader.read_variable(variable_name)
    }
    
    fn read_variable_slice(
        &self, 
        variable_name: &str, 
        indices: &[(usize, usize, usize)]
    ) -> Result<Array4<f32>, DataReaderError> {
        self.reader.read_variable_slice(variable_name, indices)
    }
    
    fn read_variables(&self, variable_names: &[&str]) -> Result<HashMap<String, Array4<f32>>, DataReaderError> {
        self.reader.read_variables(variable_names)
    }
    
    fn read_coordinates(&self) -> Result<(Array2<f32>, Array2<f32>, Array1<f32>), DataReaderError> {
        self.reader.read_coordinates()
    }
    
    fn get_global_attributes(&self) -> Result<Attributes, DataReaderError> {
        self.reader.get_global_attributes()
    }
    
    fn get_variable_attributes(&self, variable_name: &str) -> Result<Attributes, DataReaderError> {
        self.reader.get_variable_attributes(variable_name)
    }
    
    fn get_global_attribute(&self, attribute_name: &str) -> Result<AttributeValue, DataReaderError> {
        self.reader.get_global_attribute(attribute_name)
    }
    
    fn get_variable_attribute(
        &self, 
        variable_name: &str, 
        attribute_name: &str
    ) -> Result<AttributeValue, DataReaderError> {
        self.reader.get_variable_attribute(variable_name, attribute_name)
    }
    
    fn get_path(&self) -> Option<PathBuf> {
        self.reader.get_path()
    }
}

// Implement DataReader trait for Box<dyn DataReader> to enable generic usage
impl DataReader for Box<dyn DataReader> {
    fn open(&mut self, path: &Path) -> Result<(), DataReaderError> {
        self.as_mut().open(path)
    }
    
    fn is_open(&self) -> bool {
        self.as_ref().is_open()
    }
    
    fn close(&mut self) {
        self.as_mut().close()
    }
    
    fn list_variables(&self) -> Result<Vec<String>, DataReaderError> {
        self.as_ref().list_variables()
    }
    
    fn get_variable_info(&self, variable_name: &str) -> Result<VariableInfo, DataReaderError> {
        self.as_ref().get_variable_info(variable_name)
    }
    
    fn list_dimensions(&self) -> Result<Vec<String>, DataReaderError> {
        self.as_ref().list_dimensions()
    }
    
    fn get_dimension_info(&self, dimension_name: &str) -> Result<DimensionInfo, DataReaderError> {
        self.as_ref().get_dimension_info(dimension_name)
    }
    
    fn get_metadata(&self) -> Result<FileMetadata, DataReaderError> {
        self.as_ref().get_metadata()
    }
    
    fn read_variable(&self, variable_name: &str) -> Result<Array4<f32>, DataReaderError> {
        self.as_ref().read_variable(variable_name)
    }
    
    fn read_variable_slice(
        &self, 
        variable_name: &str, 
        indices: &[(usize, usize, usize)]
    ) -> Result<Array4<f32>, DataReaderError> {
        self.as_ref().read_variable_slice(variable_name, indices)
    }
    
    fn read_variables(&self, variable_names: &[&str]) -> Result<HashMap<String, Array4<f32>>, DataReaderError> {
        self.as_ref().read_variables(variable_names)
    }
    
    fn read_coordinates(&self) -> Result<(Array2<f32>, Array2<f32>, Array1<f32>), DataReaderError> {
        self.as_ref().read_coordinates()
    }
    
    fn get_global_attributes(&self) -> Result<Attributes, DataReaderError> {
        self.as_ref().get_global_attributes()
    }
    
    fn get_variable_attributes(&self, variable_name: &str) -> Result<Attributes, DataReaderError> {
        self.as_ref().get_variable_attributes(variable_name)
    }
    
    fn get_global_attribute(&self, attribute_name: &str) -> Result<AttributeValue, DataReaderError> {
        self.as_ref().get_global_attribute(attribute_name)
    }
    
    fn get_variable_attribute(
        &self, 
        variable_name: &str, 
        attribute_name: &str
    ) -> Result<AttributeValue, DataReaderError> {
        self.as_ref().get_variable_attribute(variable_name, attribute_name)
    }
    
    fn get_path(&self) -> Option<PathBuf> {
        self.as_ref().get_path()
    }
}

/// Auto-detect file format and create appropriate reader
/// 
/// This function implements the format detection strategy:
/// 1. First tries NetCDF magic bytes detection
/// 2. If that fails, checks for Zarr format indicators:
///    - `.zarr` suffix in path
///    - Directory with `.zgroup` or `.zarray` files
///    - JSON header with Zarr metadata
/// 
/// # Arguments
/// * `path` - Path to data file/dataset
/// 
/// # Returns
/// * `Result<Box<dyn DataReader>, DataReaderError>` - Appropriate reader or error
pub fn create_reader<P: AsRef<Path>>(path: P) -> Result<Box<dyn DataReader>, DataReaderError> {
    let path_ref = path.as_ref();
    
    // First, try to detect NetCDF by magic bytes
    if let Ok(true) = is_netcdf_format(path_ref) {
        return create_netcdf_reader(path);
    }
    
    // If NetCDF detection failed, check for Zarr format
    if let Ok(true) = is_zarr_format(path_ref) {
        return create_zarr_reader(path);
    }
    
    // Fallback to extension-based detection for compatibility
    if let Some(extension) = path_ref.extension() {
        match extension.to_str() {
            Some("nc") | Some("nc4") | Some("netcdf") => create_netcdf_reader(path),
            Some("zarr") => create_zarr_reader(path),
            _ => {
                Err(DataReaderError::InvalidFormat(
                    format!("Unable to detect format for file: {:?}", path_ref)
                ))
            }
        }
    } else {
        Err(DataReaderError::InvalidFormat(
            format!("Unable to detect format for file: {:?}", path_ref)
        ))
    }
}

/// Format detection utilities
/// Check if a file is in NetCDF format by examining magic bytes
/// 
/// NetCDF files start with specific magic bytes:
/// - Classic NetCDF: "CDF\001" or "CDF\002"
/// - NetCDF-4 (HDF5): "\211HDF\r\n\032\n"
/// 
/// # Arguments
/// * `path` - Path to the file to check
/// 
/// # Returns
/// * `Result<bool, DataReaderError>` - True if NetCDF format detected, false otherwise
pub fn is_netcdf_format(path: &Path) -> Result<bool, DataReaderError> {
    if !path.exists() {
        return Err(DataReaderError::FileNotFound(path.to_string_lossy().to_string()));
    }
    
    if path.is_dir() {
        return Ok(false);
    }
    
    let mut file = File::open(path)?;
    let mut buffer = [0u8; 8];
    
    // Read first 8 bytes
    match file.read_exact(&mut buffer) {
        Ok(_) => {
            // Check for classic NetCDF magic bytes
            if buffer[0..3] == [b'C', b'D', b'F'] && (buffer[3] == 1 || buffer[3] == 2) {
                return Ok(true);
            }
            
            // Check for NetCDF-4 (HDF5) magic bytes
            if buffer == [0x89, b'H', b'D', b'F', b'\r', b'\n', 0x1a, b'\n'] {
                return Ok(true);
            }
            
            Ok(false)
        },
        Err(_) => {
            // File too small or read error - not NetCDF
            Ok(false)
        }
    }
}

/// Check if a path refers to a Zarr dataset
/// 
/// This function checks for Zarr format indicators:
/// 1. `.zarr` suffix in the path name
/// 2. Directory containing `.zgroup` file (Zarr group)
/// 3. Directory containing `.zarray` file (Zarr array)
/// 4. File with JSON content containing Zarr metadata
/// 
/// # Arguments
/// * `path` - Path to check for Zarr format
/// 
/// # Returns
/// * `Result<bool, DataReaderError>` - True if Zarr format detected, false otherwise
pub fn is_zarr_format(path: &Path) -> Result<bool, DataReaderError> {
    if !path.exists() {
        return Err(DataReaderError::FileNotFound(path.to_string_lossy().to_string()));
    }
    
    // Check for .zarr suffix
    if let Some(extension) = path.extension() {
        if extension == "zarr" {
            return Ok(true);
        }
    }
    
    // Check for .zarr suffix in the file name (not just extension)
    if let Some(file_name) = path.file_name() {
        if file_name.to_string_lossy().contains(".zarr") {
            return Ok(true);
        }
    }
    
    if path.is_dir() {
        // Check for Zarr group marker
        if path.join(".zgroup").exists() {
            return Ok(true);
        }
        
        // Check for Zarr array marker
        if path.join(".zarray").exists() {
            return Ok(true);
        }
        
        // Check if any subdirectories contain Zarr markers
        if let Ok(entries) = std::fs::read_dir(path) {
            for entry in entries.flatten() {
                let entry_path = entry.path();
                if entry_path.is_dir() && (entry_path.join(".zarray").exists() || entry_path.join(".zgroup").exists()) {
                    return Ok(true);
                }
            }
        }
    } else {
        // Check if it's a JSON file with Zarr metadata
        if let Ok(true) = is_zarr_json_file(path) {
            return Ok(true);
        }
    }
    
    Ok(false)
}

/// Check if a file is a JSON file containing Zarr metadata
/// 
/// # Arguments
/// * `path` - Path to the file to check
/// 
/// # Returns
/// * `Result<bool, DataReaderError>` - True if file contains Zarr JSON metadata
fn is_zarr_json_file(path: &Path) -> Result<bool, DataReaderError> {
    let mut file = File::open(path)?;
    let mut reader = BufReader::new(&mut file);
    
    // Read a reasonable amount of the file to check for JSON content
    let mut buffer = [0u8; 1024];
    let content = match reader.read(&mut buffer) {
        Ok(bytes_read) => {
            String::from_utf8_lossy(&buffer[..bytes_read]).to_string()
        },
        Err(_) => return Ok(false),
    };
    
    // Check if content looks like JSON and contains Zarr-specific fields
    if content.trim_start().starts_with('{') {
        // Look for Zarr-specific JSON fields
        let zarr_indicators = [
            "zarr_format",
            "chunks", 
            "compressor",
            "dtype",
            "fill_value",
            "filters",
            "order",
            "shape",
            "dimension_separator"
        ];
        
        for indicator in &zarr_indicators {
            if content.contains(indicator) {
                return Ok(true);
            }
        }
    }
    
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_attribute_value_types() {
        let string_attr = AttributeValue::String("test".to_string());
        let int_attr = AttributeValue::Int(42);
        let float_attr = AttributeValue::Float(std::f32::consts::PI);
        
        match string_attr {
            AttributeValue::String(s) => assert_eq!(s, "test"),
            _ => panic!("Expected string attribute"),
        }
        
        match int_attr {
            AttributeValue::Int(i) => assert_eq!(i, 42),
            _ => panic!("Expected int attribute"),
        }
        
        match float_attr {
            AttributeValue::Float(f) => assert!((f - std::f32::consts::PI).abs() < 1e-6),
            _ => panic!("Expected float attribute"),
        }
    }
    
    #[test]
    fn test_error_types() {
        let error = DataReaderError::MissingVariable("test_var".to_string());
        assert!(error.to_string().contains("test_var"));
        
        let error = DataReaderError::FileNotFound("test.nc".to_string());
        assert!(error.to_string().contains("test.nc"));
    }
    
    // Note: Comprehensive tests for format detection would require tempfile crate
    // These tests demonstrate the basic functionality without external dependencies
    
    #[test] 
    fn test_format_detection_basic() {
        // Test nonexistent file handling
        assert!(is_netcdf_format(std::path::Path::new("/nonexistent/file.nc")).is_err());
        assert!(is_zarr_format(std::path::Path::new("/nonexistent/path")).is_err());
        
        // Test that the detection functions exist and handle errors properly
        let result = is_netcdf_format(std::path::Path::new("/tmp"));
        // Directory should return Ok(false) for NetCDF
        if result.is_ok() {
            assert!(!result.unwrap());
        }
    }
}
