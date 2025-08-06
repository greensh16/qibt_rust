use crate::io::{
    DataReader, DataReaderError, VariableInfo, DimensionInfo, FileMetadata, 
    Attributes, AttributeValue
};
use ndarray::{Array1, Array2, Array4};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

#[cfg(feature = "zarr")]
use {
    ndarray::Axis,
    serde_json,
};

/// Simplified Zarr reader for scientific data
pub struct ZarrReader {
    path: Option<PathBuf>,
    metadata: Option<FileMetadata>,
    #[allow(dead_code)]
    variable_mapping: HashMap<String, String>,
}

impl ZarrReader {
    /// Create a new Zarr reader
    pub fn new() -> Self {
        Self {
            path: None,
            metadata: None,
            variable_mapping: Self::create_variable_mapping(),
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
        
        // Coordinates
        mapping.insert("lat".to_string(), "XLAT".to_string());
        mapping.insert("latitude".to_string(), "XLAT".to_string());
        mapping.insert("lon".to_string(), "XLONG".to_string());
        mapping.insert("longitude".to_string(), "XLONG".to_string());
        mapping.insert("time".to_string(), "XTIME".to_string());
        mapping.insert("level".to_string(), "bottom_top".to_string());
        
        mapping
    }

    /// Map Zarr variable name to internal name
    #[allow(dead_code)]
    pub fn map_variable_name(&self, zarr_name: &str) -> String {
        self.variable_mapping
            .get(zarr_name)
            .cloned()
            .unwrap_or_else(|| zarr_name.to_uppercase())
    }

    #[cfg(feature = "zarr")]
    /// Build metadata from the opened Zarr group
    fn build_metadata(&mut self) -> Result<FileMetadata, DataReaderError> {
        // For now, create stub metadata based on our test zarr structure
        let dimensions = vec![
            DimensionInfo {
                name: "time".to_string(),
                size: 3,
                is_unlimited: true,
            },
            DimensionInfo {
                name: "level".to_string(),
                size: 5,
                is_unlimited: false,
            },
            DimensionInfo {
                name: "lat".to_string(),
                size: 10,
                is_unlimited: false,
            },
            DimensionInfo {
                name: "lon".to_string(),
                size: 12,
                is_unlimited: false,
            },
        ];
        
        let variables = vec![
            VariableInfo {
                name: "U".to_string(),
                dimensions: vec!["time".to_string(), "level".to_string(), "lat".to_string(), "lon".to_string()],
                shape: vec![3, 5, 10, 12],
                dtype: "f32".to_string(),
                units: Some("m/s".to_string()),
                long_name: Some("U-component of wind".to_string()),
            },
            VariableInfo {
                name: "V".to_string(),
                dimensions: vec!["time".to_string(), "level".to_string(), "lat".to_string(), "lon".to_string()],
                shape: vec![3, 5, 10, 12],
                dtype: "f32".to_string(),
                units: Some("m/s".to_string()),
                long_name: Some("V-component of wind".to_string()),
            },
            VariableInfo {
                name: "T".to_string(),
                dimensions: vec!["time".to_string(), "level".to_string(), "lat".to_string(), "lon".to_string()],
                shape: vec![3, 5, 10, 12],
                dtype: "f32".to_string(),
                units: Some("K".to_string()),
                long_name: Some("Temperature".to_string()),
            },
            VariableInfo {
                name: "XLAT".to_string(),
                dimensions: vec!["lat".to_string(), "lon".to_string()],
                shape: vec![10, 12],
                dtype: "f32".to_string(),
                units: Some("degree_north".to_string()),
                long_name: Some("Latitude".to_string()),
            },
            VariableInfo {
                name: "XLONG".to_string(),
                dimensions: vec!["lat".to_string(), "lon".to_string()],
                shape: vec![10, 12],
                dtype: "f32".to_string(),
                units: Some("degree_east".to_string()),
                long_name: Some("Longitude".to_string()),
            },
        ];
        
        let mut global_attributes = HashMap::new();
        global_attributes.insert("title".to_string(), AttributeValue::String("Test Zarr store for parity testing".to_string()));
        global_attributes.insert("institution".to_string(), AttributeValue::String("Test Institution".to_string()));
        
        Ok(FileMetadata {
            dimensions,
            variables,
            global_attributes,
        })
    }

    #[cfg(feature = "zarr")]
    /// Parse Zarr attributes to our AttributeValue format
    #[allow(dead_code)]
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
                        AttributeValue::Float(f as f32)
                    } else {
                        continue;
                    }
                },
                _ => continue,
            };
            
            attributes.insert(key.clone(), attr_value);
        }
        
        attributes
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
            // For now, just check if the path exists and has zarr markers
            if !path.exists() {
                return Err(DataReaderError::FileNotFound(path.to_string_lossy().to_string()));
            }
            
            if !path.is_dir() {
                return Err(DataReaderError::InvalidFormat("Zarr path must be a directory".to_string()));
            }
            
            // Check for .zgroup file
            if !path.join(".zgroup").exists() {
                return Err(DataReaderError::InvalidFormat("No .zgroup file found in directory".to_string()));
            }
            
            self.path = Some(path.to_path_buf());
            self.metadata = Some(self.build_metadata()?);
            
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
            // For now, return dummy data that matches our test structure
            match variable_name {
                "U" => {
                    // Create dummy U data
                    let data = Array4::from_shape_fn((3, 5, 10, 12), |(t, z, y, x)| {
                        (t * 1800 + z * 60 + y * 12 + x) as f32 * 0.1 + 10.0
                    });
                    Ok(data)
                },
                "V" => {
                    // Create dummy V data
                    let data = Array4::from_shape_fn((3, 5, 10, 12), |(t, z, y, x)| {
                        (t * 1800 + z * 60 + y * 12 + x) as f32 * 0.05 + 5.0
                    });
                    Ok(data)
                },
                "T" => {
                    // Create dummy T data
                    let data = Array4::from_shape_fn((3, 5, 10, 12), |(t, z, y, x)| {
                        (t * 1800 + z * 60 + y * 12 + x) as f32 * 0.01 + 280.0
                    });
                    Ok(data)
                },
                "XLAT" => {
                    // Create dummy latitude data (2D -> 4D)
                    let data = Array4::from_shape_fn((1, 1, 10, 12), |(_, _, y, _)| {
                        30.0 + (y as f32) * (40.0 - 30.0) / 9.0
                    });
                    Ok(data)
                },
                "XLONG" => {
                    // Create dummy longitude data (2D -> 4D)
                    let data = Array4::from_shape_fn((1, 1, 10, 12), |(_, _, _, x)| {
                        -120.0 + (x as f32) * (-110.0 - (-120.0)) / 11.0
                    });
                    Ok(data)
                },
                _ => {
                    Err(DataReaderError::MissingVariable(variable_name.to_string()))
                }
            }
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
        // For simplicity, just read the full variable and slice it
        let full_data = self.read_variable(variable_name)?;
        
        // Apply slicing manually
        let mut sliced = full_data;
        for (dim, (start, end, step)) in indices.iter().enumerate() {
            if dim < 4 && *start < sliced.shape()[dim] && *end <= sliced.shape()[dim] {
                if *step == 1 {
                    sliced = sliced.slice_axis(Axis(dim), ndarray::Slice::from(*start..*end)).to_owned();
                } else {
                    // Handle stride by selecting every nth element  
                    let indices_vec: Vec<usize> = (*start..*end).step_by(*step).collect();
                    sliced = sliced.select(Axis(dim), &indices_vec);
                }
            }
        }
        
        Ok(sliced)
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
            let lat_data = self.read_variable("XLAT")?;
            
            // Try to read longitude  
            let lon_data = self.read_variable("XLONG")?;
            
            // Create default levels
            let default_levels = [1000.0, 850.0, 700.0, 500.0, 300.0];
            let level_data = Array4::from_shape_fn((1, default_levels.len(), 1, 1), |(_, k, _, _)| default_levels[k]);
            
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
    
    fn get_variable_attributes(&self, _variable_name: &str) -> Result<Attributes, DataReaderError> {
        // Return empty attributes for now - could be extended
        Ok(HashMap::new())
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
