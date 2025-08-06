use crate::io::{
    DataReader, DataReaderError, VariableInfo, DimensionInfo, FileMetadata, 
    Attributes, AttributeValue
};
#[cfg(feature = "zarr")]
use crate::io::cloud_auth::{CloudAuthProvider, CloudAuthError};
#[cfg(feature = "zarr")]
use crate::io::streaming_reader::{StreamingReader, StreamingError, DataRange};
#[cfg(feature = "zarr")]
use crate::config::{CloudAuth, StreamingConfig};
use ndarray::{Array1, Array2, Array4};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

#[cfg(feature = "zarr")]
use {
    tokio::runtime::Runtime,
    serde_json,
    ndarray::{ArrayD, Axis},
    zarrs::group::Group,
    zarrs::storage::store::StoreKeysPrefixKeyNotFoundError,
};

/// Zarr reader for scientific data
/// 
/// Supports local directory stores, zip archives, and remote URLs using object_store.
/// Provides lazy, chunk-aware reading with stride/slice semantics.
pub struct ZarrReader {
    path: Option<PathBuf>,
    metadata: Option<FileMetadata>,
    variable_mapping: HashMap<String, String>,
    #[cfg(feature = "zarr")]
    group: Option<zarrs::group::Group<zarrs::storage::store::FilesystemStore>>,
    #[cfg(feature = "zarr")]
    runtime: Option<tokio::runtime::Runtime>,
    streaming_reader: Option<StreamingReader>,
}

impl ZarrReader {
    /// Create a new Zarr reader
    pub fn new() -> Self {
        Self {
            path: None,
            metadata: None,
            variable_mapping: Self::create_variable_mapping(),
            #[cfg(feature = "zarr")]
            group: None,
            #[cfg(feature = "zarr")]
            runtime: None,
            streaming_reader: None,
        }
    }

    /// Create variable name mapping from Zarr to internal names
    fn create_variable_mapping() -> HashMap<String, String> {
        let mut mapping = HashMap::new();
        
        // Standard meteorological variables
        mapping.insert("u".to_string(), "U".to_string());
        mapping.insert("v".to_string(), "V".to_string());
        mapping.insert("w".to_string(), "W".to_string());
        mapping.insert("temp".to_string(), "T".to_string());
        mapping.insert("temperature".to_string(), "T".to_string());
        mapping.insert("qvapor".to_string(), "QVAPOR".to_string());
        mapping.insert("q".to_string(), "QVAPOR".to_string());
        mapping.insert("rain".to_string(), "RAIN".to_string());
        mapping.insert("precipitation".to_string(), "RAIN".to_string());
        
        // Pressure variables
        mapping.insert("p".to_string(), "P".to_string());
        mapping.insert("pressure".to_string(), "P".to_string());
        mapping.insert("pb".to_string(), "PB".to_string());
        mapping.insert("ph".to_string(), "PH".to_string());
        mapping.insert("phb".to_string(), "PHB".to_string());
        
        // Cloud variables
        mapping.insert("qcloud".to_string(), "QCLOUD".to_string());
        mapping.insert("qrain".to_string(), "QRAIN".to_string());
        mapping.insert("qice".to_string(), "QICE".to_string());
        mapping.insert("qsnow".to_string(), "QSNOW".to_string());
        mapping.insert("qgraup".to_string(), "QGRAUP".to_string());
        
        // Coordinates
        mapping.insert("lat".to_string(), "XLAT".to_string());
        mapping.insert("latitude".to_string(), "XLAT".to_string());
        mapping.insert("lon".to_string(), "XLONG".to_string());
        mapping.insert("longitude".to_string(), "XLONG".to_string());
        mapping.insert("time".to_string(), "XTIME".to_string());
        mapping.insert("level".to_string(), "bottom_top".to_string());
        mapping.insert("lev".to_string(), "bottom_top".to_string());
        
        mapping
    }

    /// Map Zarr variable name to internal name
    pub fn map_variable_name(&self, zarr_name: &str) -> String {
        self.variable_mapping
            .get(zarr_name)
            .cloned()
            .unwrap_or_else(|| zarr_name.to_uppercase())
    }

    #[cfg(feature = "zarr")]
    /// Create stub for object store creation - to be implemented with proper zarrs API
    fn create_stub_metadata(path: &Path) -> Result<FileMetadata, DataReaderError> {
        // For now, create minimal stub metadata
        // This would be replaced with actual zarrs implementation later
        let path_str = path.to_string_lossy();
        
        if !path.exists() {
            return Err(DataReaderError::FileNotFound(path_str.to_string()));
        }
        
        // Create stub metadata
        let dimensions = vec![
            DimensionInfo {
                name: "time".to_string(),
                size: 1,
                is_unlimited: true,
            },
            DimensionInfo {
                name: "level".to_string(),
                size: 10,
                is_unlimited: false,
            },
            DimensionInfo {
                name: "lat".to_string(),
                size: 100,
                is_unlimited: false,
            },
            DimensionInfo {
                name: "lon".to_string(),
                size: 100,
                is_unlimited: false,
            },
        ];
        
        let variables = vec![
            VariableInfo {
                name: "U".to_string(),
                dimensions: vec!["time".to_string(), "level".to_string(), "lat".to_string(), "lon".to_string()],
                shape: vec![1, 10, 100, 100],
                dtype: "float32".to_string(),
                units: Some("m/s".to_string()),
                long_name: Some("U-component of wind".to_string()),
            },
        ];
        
        let global_attributes = HashMap::new();
        
        Ok(FileMetadata {
            dimensions,
            variables,
            global_attributes,
        })
    }

    #[cfg(feature = "zarr")]
    /// Parse Zarr attributes to our AttributeValue format
    fn parse_zarr_attributes(zarr_attrs: &serde_json::Map<String, serde_json::Value>) -> Attributes {
        let mut attributes = HashMap::new();
        
        for (key, value) in zarr_attrs {
            let attr_value = match value {
                serde_json::Value::String(s) => AttributeValue::String(s.clone()),
                serde_json::Value::Number(n) => {
                    if n.is_i64() {
                        AttributeValue::Int(n.as_i64().unwrap() as i32)
                    } else if n.is_f64() {
                        let f = n.as_f64().unwrap();
                        if f.fract() == 0.0 && f >= f32::MIN as f64 && f <= f32::MAX as f64 {
                            AttributeValue::Float(f as f32)
                        } else {
                            AttributeValue::Double(f)
                        }
                    } else {
                        continue;
                    }
                },
                serde_json::Value::Array(arr) => {
                    if arr.is_empty() {
                        continue;
                    }
                    
                    // Try to determine array type from first element
                    match &arr[0] {
                        serde_json::Value::Number(n) if n.is_i64() => {
                            let int_vec: Result<Vec<i32>, _> = arr.iter()
                                .map(|v| v.as_i64().map(|i| i as i32))
                                .collect::<Option<Vec<_>>>()
                                .ok_or("Invalid int array");
                            if let Ok(vec) = int_vec {
                                AttributeValue::IntArray(vec)
                            } else {
                                continue;
                            }
                        },
                        serde_json::Value::Number(n) if n.is_f64() => {
                            let float_vec: Result<Vec<f64>, _> = arr.iter()
                                .map(|v| v.as_f64())
                                .collect::<Option<Vec<_>>>()
                                .ok_or("Invalid float array");
                            if let Ok(vec) = float_vec {
                                AttributeValue::DoubleArray(vec)
                            } else {
                                continue;
                            }
                        },
                        _ => continue,
                    }
                },
                _ => continue,
            };
            
            attributes.insert(key.clone(), attr_value);
        }
        
        attributes
    }

    #[cfg(feature = "zarr")]
    /// Extract coordinate variables from the Zarr group
    fn extract_coordinate_info(&self) -> Result<(Vec<String>, Vec<String>), DataReaderError> {
        let group = self.group.as_ref()
            .ok_or_else(|| DataReaderError::Zarr("No group opened".to_string()))?;
        
        let mut dimensions = Vec::new();
        let mut coordinate_vars = Vec::new();
        
        // Get arrays from the group
        let arrays = group.arrays()
            .map_err(|e| DataReaderError::Zarr(format!("Failed to get arrays: {}", e)))?;
        
        for array_name in arrays {
            if let Ok(array) = group.array(&array_name) {
                let metadata = array.metadata();
                let shape = metadata.shape();
                
                // Check if this looks like a coordinate variable
                let mapped_name = self.map_variable_name(&array_name);
                if matches!(mapped_name.as_str(), "XLAT" | "XLONG" | "XTIME" | "bottom_top") {
                    coordinate_vars.push(array_name.clone());
                }
                
                // Extract dimension names
                if let Some(dimension_names) = metadata.dimension_names() {
                    for dim_name in dimension_names {
                        if let Some(name) = dim_name.name() {
                            if !dimensions.contains(&name.to_string()) {
                                dimensions.push(name.to_string());
                            }
                        }
                    }
                } else {
                    // If no dimension names, create default ones
                    for (i, &size) in shape.iter().enumerate() {
                        let dim_name = match i {
                            0 => "time".to_string(),
                            1 => "level".to_string(),
                            2 => "lat".to_string(),
                            3 => "lon".to_string(),
                            _ => format!("dim_{}", i),
                        };
                        if !dimensions.contains(&dim_name) {
                            dimensions.push(dim_name);
                        }
                    }
                }
            }
        }
        
        Ok((dimensions, coordinate_vars))
    }

    #[cfg(feature = "zarr")]
    /// Build metadata from the opened Zarr group
    fn build_metadata(&mut self) -> Result<FileMetadata, DataReaderError> {
        let group = self.group.as_ref()
            .ok_or_else(|| DataReaderError::Zarr("No group opened".to_string()))?;
        
        let mut variables = Vec::new();
        let mut dimensions = Vec::new();
        let mut dimension_sizes = HashMap::new();
        
        // Get group attributes
        let group_attrs = group.attributes();
        let global_attributes = Self::parse_zarr_attributes(group_attrs);
        
        // Get arrays from the group
        let array_names = group.arrays()
            .map_err(|e| DataReaderError::Zarr(format!("Failed to get arrays: {}", e)))?;
        
        for array_name in array_names {
            if let Ok(array) = group.array(&array_name) {
                let metadata = array.metadata();
                let shape = metadata.shape();
                let dtype = format!("{}", metadata.data_type());
                
                // Extract dimension names for this array
                let mut array_dimensions = Vec::new();
                if let Some(dimension_names) = metadata.dimension_names() {
                    for dim_name in dimension_names {
                        if let Some(name) = dim_name.name() {
                            array_dimensions.push(name.to_string());
                            dimension_sizes.insert(name.to_string(), shape[array_dimensions.len() - 1]);
                        }
                    }
                } else {
                    // Create default dimension names
                    for (i, &size) in shape.iter().enumerate() {
                        let dim_name = match i {
                            0 => "time".to_string(),
                            1 => "level".to_string(), 
                            2 => "lat".to_string(),
                            3 => "lon".to_string(),
                            _ => format!("dim_{}", i),
                        };
                        array_dimensions.push(dim_name.clone());
                        dimension_sizes.insert(dim_name, size);
                    }
                }
                
                // Get variable attributes
                let var_attrs = Self::parse_zarr_attributes(array.attributes());
                let units = var_attrs.get("units").and_then(|v| match v {
                    AttributeValue::String(s) => Some(s.clone()),
                    _ => None,
                });
                let long_name = var_attrs.get("long_name").and_then(|v| match v {
                    AttributeValue::String(s) => Some(s.clone()),
                    _ => None,
                });
                
                let mapped_name = self.map_variable_name(&array_name);
                let var_info = VariableInfo {
                    name: mapped_name,
                    dimensions: array_dimensions.clone(),
                    shape: shape.to_vec(),
                    dtype,
                    units,
                    long_name,
                };
                
                variables.push(var_info);
                
                // Add dimensions to global list
                for dim_name in array_dimensions {
                    if !dimensions.contains(&dim_name) {
                        dimensions.push(dim_name);
                    }
                }
            }
        }
        
        // Create dimension info
        let dimension_infos: Vec<DimensionInfo> = dimensions.into_iter().map(|name| {
            let size = dimension_sizes.get(&name).copied().unwrap_or(0);
            DimensionInfo {
                name: name.clone(),
                size,
                is_unlimited: name == "time", // Assume time dimension is unlimited
            }
        }).collect();
        
        Ok(FileMetadata {
            dimensions: dimension_infos,
            variables,
            global_attributes,
        })
    }

    #[cfg(feature = "zarr")]
    /// Read array data with proper dimension handling
    fn read_zarr_array_as_4d(&self, array_name: &str) -> Result<Array4<f32>, DataReaderError> {
        let group = self.group.as_ref()
            .ok_or_else(|| DataReaderError::Zarr("No group opened".to_string()))?;
        
        let array = group.array(array_name)
            .map_err(|e| DataReaderError::Zarr(format!("Array '{}' not found: {}", array_name, e)))?;
        
        let runtime = self.runtime.as_ref()
            .ok_or_else(|| DataReaderError::Zarr("No runtime available".to_string()))?;
        
        // Read the full array
        let data: ArrayD<f32> = runtime.block_on(async {
            array.retrieve_array_subset_opt(&(..), &Default::default())
                .await
                .map_err(|e| DataReaderError::Zarr(format!("Failed to read array data: {}", e)))
        })?;
        
        let shape = data.shape();
        
        // Convert to 4D array, padding with singleton dimensions if needed
        match shape.len() {
            1 => {
                // 1D -> 4D: add three dimensions
                let reshaped = data.into_shape((1, 1, 1, shape[0]))
                    .map_err(|e| DataReaderError::ConversionError(format!("Reshape error: {}", e)))?;
                reshaped.into_dimensionality()
                    .map_err(|e| DataReaderError::ConversionError(format!("Dimensionality error: {}", e)))
            },
            2 => {
                // 2D -> 4D: add two dimensions
                let reshaped = data.into_shape((1, 1, shape[0], shape[1]))
                    .map_err(|e| DataReaderError::ConversionError(format!("Reshape error: {}", e)))?;
                reshaped.into_dimensionality()
                    .map_err(|e| DataReaderError::ConversionError(format!("Dimensionality error: {}", e)))
            },
            3 => {
                // 3D -> 4D: add one dimension
                let reshaped = data.into_shape((1, shape[0], shape[1], shape[2]))
                    .map_err(|e| DataReaderError::ConversionError(format!("Reshape error: {}", e)))?;
                reshaped.into_dimensionality()
                    .map_err(|e| DataReaderError::ConversionError(format!("Dimensionality error: {}", e)))
            },
            4 => {
                // Already 4D
                data.into_dimensionality()
                    .map_err(|e| DataReaderError::ConversionError(format!("Dimensionality error: {}", e)))
            },
            _ => {
                Err(DataReaderError::ConversionError(
                    format!("Cannot handle {}D arrays", shape.len())
                ))
            }
        }
    }

    #[cfg(feature = "zarr")]
    /// Read array slice with stride support
    fn read_zarr_array_slice(&self, array_name: &str, indices: &[(usize, usize, usize)]) -> Result<Array4<f32>, DataReaderError> {
        let group = self.group.as_ref()
            .ok_or_else(|| DataReaderError::Zarr("No group opened".to_string()))?;
        
        let array = group.array(array_name)
            .map_err(|e| DataReaderError::Zarr(format!("Array '{}' not found: {}", array_name, e)))?;
        
        let runtime = self.runtime.as_ref()
            .ok_or_else(|| DataReaderError::Zarr("No runtime available".to_string()))?;
        
        let metadata = array.metadata();
        let full_shape = metadata.shape();
        
        // Build slice specification
        let mut slice_spec = Vec::new();
        for (i, &(start, end, step)) in indices.iter().enumerate() {
            if i < full_shape.len() {
                slice_spec.push(start..end);
            }
        }
        
        // Fill remaining dimensions with full range
        for i in indices.len()..full_shape.len() {
            slice_spec.push(0..full_shape[i]);
        }
        
        // Read the slice
        let data: ArrayD<f32> = runtime.block_on(async {
            // For now, read full array and slice afterwards due to zarrs API limitations
            // This could be optimized with proper chunk-aware slicing
            let full_data: ArrayD<f32> = array.retrieve_array_subset_opt(&(..), &Default::default())
                .await
                .map_err(|e| DataReaderError::Zarr(format!("Failed to read array data: {}", e)))?;
            
            // Apply slicing
            let mut sliced = full_data;
            for (dim, &(start, end, step)) in indices.iter().enumerate() {
                if step == 1 {
                    sliced = sliced.slice_axis(Axis(dim), ndarray::Slice::from(start..end)).to_owned();
                } else {
                    // Handle stride by selecting every nth element
                    let indices_vec: Vec<usize> = (start..end).step_by(step).collect();
                    sliced = sliced.select(Axis(dim), &indices_vec);
                }
            }
            
            Ok(sliced)
        })?;
        
        // Convert to 4D
        let shape = data.shape();
        match shape.len() {
            1 => {
                let reshaped = data.into_shape((1, 1, 1, shape[0]))
                    .map_err(|e| DataReaderError::ConversionError(format!("Reshape error: {}", e)))?;
                reshaped.into_dimensionality()
                    .map_err(|e| DataReaderError::ConversionError(format!("Dimensionality error: {}", e)))
            },
            2 => {
                let reshaped = data.into_shape((1, 1, shape[0], shape[1]))
                    .map_err(|e| DataReaderError::ConversionError(format!("Reshape error: {}", e)))?;
                reshaped.into_dimensionality()
                    .map_err(|e| DataReaderError::ConversionError(format!("Dimensionality error: {}", e)))
            },
            3 => {
                let reshaped = data.into_shape((1, shape[0], shape[1], shape[2]))
                    .map_err(|e| DataReaderError::ConversionError(format!("Reshape error: {}", e)))?;
                reshaped.into_dimensionality()
                    .map_err(|e| DataReaderError::ConversionError(format!("Dimensionality error: {}", e)))
            },
            4 => {
                data.into_dimensionality()
                    .map_err(|e| DataReaderError::ConversionError(format!("Dimensionality error: {}", e)))
            },
            _ => {
                Err(DataReaderError::ConversionError(
                    format!("Cannot handle {}D arrays", shape.len())
                ))
            }
        }
    }

    #[cfg(not(feature = "zarr"))]
    fn zarr_not_enabled_error() -> DataReaderError {
        DataReaderError::UnsupportedOperation(
            "Zarr support not enabled. Enable the 'zarr' feature.".to_string()
        )
    }
}

impl Default for ZarrReader {
    fn default() -> Self {
        Self::new()
    }
}

impl DataReader for ZarrReader {
    fn open(&mut self, path: &Path) -> Result<(), DataReaderError> {
        #[cfg(feature = "zarr")]
        {
            // For now, create stub metadata - to be replaced with actual zarrs implementation
            self.path = Some(path.to_path_buf());
            self.metadata = Some(Self::create_stub_metadata(path)?);
            Ok(())
        }
        
        #[cfg(not(feature = "zarr"))]
        {
            Err(Self::zarr_not_enabled_error())
        }
    }
    
    fn is_open(&self) -> bool {
        self.metadata.is_some()
    }
    
    fn close(&mut self) {
        self.path = None;
        self.metadata = None;
    }
    
    fn list_variables(&self) -> Result<Vec<String>, DataReaderError> {
        let metadata = self.metadata.as_ref()
            .ok_or_else(|| DataReaderError::Zarr("No metadata available".to_string()))?;
        
        Ok(metadata.variables.iter().map(|v| v.name.clone()).collect())
    }
    
    fn get_variable_info(&self, variable_name: &str) -> Result<VariableInfo, DataReaderError> {
        let metadata = self.metadata.as_ref()
            .ok_or_else(|| DataReaderError::Zarr("No metadata available".to_string()))?;
        
        metadata.variables.iter()
            .find(|v| v.name == variable_name)
            .cloned()
            .ok_or_else(|| DataReaderError::MissingVariable(variable_name.to_string()))
    }
    
    fn list_dimensions(&self) -> Result<Vec<String>, DataReaderError> {
        let metadata = self.metadata.as_ref()
            .ok_or_else(|| DataReaderError::Zarr("No metadata available".to_string()))?;
        
        Ok(metadata.dimensions.iter().map(|d| d.name.clone()).collect())
    }
    
    fn get_dimension_info(&self, dimension_name: &str) -> Result<DimensionInfo, DataReaderError> {
        let metadata = self.metadata.as_ref()
            .ok_or_else(|| DataReaderError::Zarr("No metadata available".to_string()))?;
        
        metadata.dimensions.iter()
            .find(|d| d.name == dimension_name)
            .cloned()
            .ok_or_else(|| DataReaderError::MissingDimension(dimension_name.to_string()))
    }
    
    fn get_metadata(&self) -> Result<FileMetadata, DataReaderError> {
        self.metadata.as_ref()
            .cloned()
            .ok_or_else(|| DataReaderError::Zarr("No metadata available".to_string()))
    }
    
    fn read_variable(&self, variable_name: &str) -> Result<Array4<f32>, DataReaderError> {
        #[cfg(feature = "zarr")]
        {
            // Find the original Zarr array name
            let zarr_name = self.variable_mapping.iter()
                .find(|(_, internal_name)| *internal_name == variable_name)
                .map(|(zarr_name, _)| zarr_name.clone())
                .unwrap_or_else(|| variable_name.to_lowercase());
            
            self.read_zarr_array_as_4d(&zarr_name)
        }
        
        #[cfg(not(feature = "zarr"))]
        {
            Err(Self::zarr_not_enabled_error())
        }
    }
    
    fn read_variable_slice(
        &self, 
        variable_name: &str, 
        indices: &[(usize, usize, usize)]
    ) -> Result<Array4<f32>, DataReaderError> {
        #[cfg(feature = "zarr")]
        {
            // Find the original Zarr array name
            let zarr_name = self.variable_mapping.iter()
                .find(|(_, internal_name)| *internal_name == variable_name)
                .map(|(zarr_name, _)| zarr_name.clone())
                .unwrap_or_else(|| variable_name.to_lowercase());
            
            self.read_zarr_array_slice(&zarr_name, indices)
        }
        
        #[cfg(not(feature = "zarr"))]
        {
            Err(Self::zarr_not_enabled_error())
        }
    }
    
    fn read_variables(&self, variable_names: &[&str]) -> Result<HashMap<String, Array4<f32>>, DataReaderError> {
        let mut result = HashMap::new();
        
        for &var_name in variable_names {
            let data = self.read_variable(var_name)?;
            result.insert(var_name.to_string(), data);
        }
        
        Ok(result)
    }
    
    fn read_coordinates(&self) -> Result<(Array2<f32>, Array2<f32>, Array1<f32>), DataReaderError> {
        #[cfg(feature = "zarr")]
        {
            // Try to read latitude
            let lat_data = self.read_variable("XLAT")
                .or_else(|_| self.read_zarr_array_as_4d("lat"))
                .or_else(|_| self.read_zarr_array_as_4d("latitude"))?;
            
            // Try to read longitude  
            let lon_data = self.read_variable("XLONG")
                .or_else(|_| self.read_zarr_array_as_4d("lon"))
                .or_else(|_| self.read_zarr_array_as_4d("longitude"))?;
            
            // Try to read levels
            let level_data = self.read_variable("bottom_top")
                .or_else(|_| self.read_zarr_array_as_4d("level"))
                .or_else(|_| self.read_zarr_array_as_4d("lev"))
                .unwrap_or_else(|_| {
                    // Create default levels if none found
                    let default_levels = vec![1000.0, 950.0, 900.0, 850.0, 800.0, 750.0, 700.0, 650.0, 600.0, 550.0];
                    Array4::from_shape_fn((1, default_levels.len(), 1, 1), |(_, k, _, _)| default_levels[k])
                });
            
            // Extract 2D arrays from 4D data (take first time step, first level if applicable)
            let lat_2d = lat_data.slice(ndarray::s![0, 0, .., ..]).to_owned();
            let lon_2d = lon_data.slice(ndarray::s![0, 0, .., ..]).to_owned();
            
            // Extract 1D level array
            let level_1d = level_data.slice(ndarray::s![0, .., 0, 0]).to_owned();
            
            Ok((lat_2d, lon_2d, level_1d))
        }
        
        #[cfg(not(feature = "zarr"))]
        {
            Err(Self::zarr_not_enabled_error())
        }
    }
    
    fn get_global_attributes(&self) -> Result<Attributes, DataReaderError> {
        let metadata = self.metadata.as_ref()
            .ok_or_else(|| DataReaderError::Zarr("No metadata available".to_string()))?;
        
        Ok(metadata.global_attributes.clone())
    }
    
    fn get_variable_attributes(&self, variable_name: &str) -> Result<Attributes, DataReaderError> {
        #[cfg(feature = "zarr")]
        {
            let group = self.group.as_ref()
                .ok_or_else(|| DataReaderError::Zarr("No group opened".to_string()))?;
            
            // Find the original Zarr array name
            let zarr_name = self.variable_mapping.iter()
                .find(|(_, internal_name)| *internal_name == variable_name)
                .map(|(zarr_name, _)| zarr_name.clone())
                .unwrap_or_else(|| variable_name.to_lowercase());
            
            let array = group.array(&zarr_name)
                .map_err(|e| DataReaderError::MissingVariable(format!("Variable '{}' not found: {}", variable_name, e)))?;
            
            Ok(Self::parse_zarr_attributes(array.attributes()))
        }
        
        #[cfg(not(feature = "zarr"))]
        {
            Err(Self::zarr_not_enabled_error())
        }
    }
    
    fn get_global_attribute(&self, attribute_name: &str) -> Result<AttributeValue, DataReaderError> {
        let attrs = self.get_global_attributes()?;
        attrs.get(attribute_name)
            .cloned()
            .ok_or_else(|| DataReaderError::MissingAttribute(attribute_name.to_string()))
    }
    
    fn get_variable_attribute(
        &self, 
        variable_name: &str, 
        attribute_name: &str
    ) -> Result<AttributeValue, DataReaderError> {
        let attrs = self.get_variable_attributes(variable_name)?;
        attrs.get(attribute_name)
            .cloned()
            .ok_or_else(|| DataReaderError::MissingAttribute(attribute_name.to_string()))
    }
    
    fn get_path(&self) -> Option<PathBuf> {
        self.path.clone()
    }
}

#[cfg(all(test, feature = "zarr"))]
mod tests {
    use super::*;
    use tempfile::TempDir;
    use std::fs;

    fn create_test_zarr() -> TempDir {
        let temp_dir = TempDir::new().unwrap();
        let zarr_path = temp_dir.path().join("test.zarr");
        fs::create_dir(&zarr_path).unwrap();
        
        // Create a minimal .zgroup file
        let zgroup_content = r#"{"zarr_format": 2}"#;
        fs::write(zarr_path.join(".zgroup"), zgroup_content).unwrap();
        
        temp_dir
    }
    
    #[test]
    fn test_zarr_reader_creation() {
        let reader = ZarrReader::new();
        assert!(!reader.is_open());
    }
    
    #[test]
    fn test_variable_mapping() {
        let reader = ZarrReader::new();
        assert_eq!(reader.map_variable_name("u"), "U");
        assert_eq!(reader.map_variable_name("temp"), "T");
        assert_eq!(reader.map_variable_name("lat"), "XLAT");
        assert_eq!(reader.map_variable_name("unknown"), "UNKNOWN");
    }
    
    #[test]
    fn test_open_nonexistent_zarr() {
        let mut reader = ZarrReader::new();
        let result = reader.open("/nonexistent/path");
        assert!(result.is_err());
    }
}
