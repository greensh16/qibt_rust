use crate::math::interpolate::{FieldDataAccess, ArrayFieldData};
use crate::data_io::{NetCDFReader, MeteoFieldArray};
use crate::io::DataReader;
use ndarray::Array4;
use std::collections::HashMap;

/// Generic wrapper that provides FieldDataAccess interface for any data reader
/// This eliminates NetCDF-specific assumptions in interpolation and trajectory code
pub struct GenericFieldAccessor {
    /// Field data in [j, i, k, t] layout
    pub data: Array4<f32>,
    /// Coordinate arrays
    pub longitudes: Vec<f64>,
    pub latitudes: Vec<f64>,
    pub levels: Vec<f64>,
    pub times: Vec<f64>,
    /// Field metadata
    pub field_name: String,
    pub units: String,
    pub missing_value: f32,
}

impl GenericFieldAccessor {
    /// Create accessor from NetCDF reader (legacy path)
    pub fn from_netcdf_reader(
        reader: &NetCDFReader,
        field_name: &str,
        coordinates: (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), // (lon, lat, lev, time)
    ) -> Result<Self, String> {
        // Read field from NetCDF reader
        let field = reader.read_field(field_name)?;
        
        // Convert legacy data to Array4
        // This is a placeholder - in real implementation, this would read actual NetCDF data
        let dummy_data = Array4::<f32>::zeros((coordinates.1.len(), coordinates.0.len(), coordinates.2.len(), coordinates.3.len()));
        
        Ok(GenericFieldAccessor {
            data: dummy_data,
            longitudes: coordinates.0,
            latitudes: coordinates.1,
            levels: coordinates.2,
            times: coordinates.3,
            field_name: field.name,
            units: field.units,
            missing_value: field.missing_value as f32,
        })
    }
    
    /// Create accessor from modern Array4 data (unified path)
    pub fn from_array4(
        data: Array4<f32>,
        coordinates: (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), // (lon, lat, lev, time)
        field_name: String,
        units: String,
        missing_value: f32,
    ) -> Self {
        GenericFieldAccessor {
            data,
            longitudes: coordinates.0,
            latitudes: coordinates.1, 
            levels: coordinates.2,
            times: coordinates.3,
            field_name,
            units,
            missing_value,
        }
    }
    
    /// Create accessor from any DataReader implementation
    pub fn from_data_reader<T: DataReader>(
        reader: &T,
        field_name: &str,
    ) -> Result<Self, String> {
        // Read the variable using the generic DataReader interface
        let data = reader.read_variable(field_name)
            .map_err(|e| format!("Failed to read variable {}: {}", field_name, e))?;
            
        // Read coordinates - this is a placeholder implementation
        // In practice, coordinates would be read from the DataReader
        let (lat_grid, lon_grid, levels) = reader.read_coordinates()
            .map_err(|e| format!("Failed to read coordinates: {}", e))?;
            
        // Extract coordinate arrays from 2D grids (simplified - assumes regular grid)
        let latitudes: Vec<f64> = lat_grid.column(0).iter().map(|&x| x as f64).collect();
        let longitudes: Vec<f64> = lon_grid.row(0).iter().map(|&x| x as f64).collect();
        let levels_vec: Vec<f64> = levels.iter().map(|&x| x as f64).collect();
        
        // Create dummy time coordinates - in real implementation this would come from reader
        let times = vec![0.0, 1.0, 2.0, 3.0]; // hours since reference
        
        // Get field attributes
        let var_attrs = reader.get_variable_attributes(field_name)
            .map_err(|e| format!("Failed to get attributes: {}", e))?;
            
        let units = var_attrs.get("units")
            .map(|attr| format!("{:?}", attr))
            .unwrap_or_else(|| "unknown".to_string());
            
        let missing_value = var_attrs.get("_FillValue")
            .and_then(|attr| match attr {
                crate::io::AttributeValue::Float(f) => Some(*f),
                crate::io::AttributeValue::Double(d) => Some(*d as f32),
                _ => None,
            })
            .unwrap_or(-9999.0f32);
        
        Ok(GenericFieldAccessor {
            data,
            longitudes,
            latitudes,
            levels: levels_vec,
            times,
            field_name: field_name.to_string(),
            units,
            missing_value,
        })
    }
    
    /// Create multiple field accessors efficiently
    pub fn from_data_reader_multi<T: DataReader>(
        reader: &T,
        field_names: &[&str],
    ) -> Result<HashMap<String, GenericFieldAccessor>, String> {
        let mut accessors = HashMap::new();
        
        // Read coordinates once for all fields
        let (lat_grid, lon_grid, levels) = reader.read_coordinates()
            .map_err(|e| format!("Failed to read coordinates: {}", e))?;
            
        let latitudes: Vec<f64> = lat_grid.column(0).iter().map(|&x| x as f64).collect();
        let longitudes: Vec<f64> = lon_grid.row(0).iter().map(|&x| x as f64).collect(); 
        let levels_vec: Vec<f64> = levels.iter().map(|&x| x as f64).collect();
        let times = vec![0.0, 1.0, 2.0, 3.0]; // Placeholder
        
        // Read all variables at once for efficiency
        let field_data = reader.read_variables(field_names)
            .map_err(|e| format!("Failed to read variables: {}", e))?;
            
        for field_name in field_names {
            if let Some(data) = field_data.get(*field_name) {
                let var_attrs = reader.get_variable_attributes(field_name)
                    .map_err(|e| format!("Failed to get attributes for {}: {}", field_name, e))?;
                    
                let units = var_attrs.get("units")
                    .map(|attr| format!("{:?}", attr))
                    .unwrap_or_else(|| "unknown".to_string());
                    
                let missing_value = var_attrs.get("_FillValue")
                    .and_then(|attr| match attr {
                        crate::io::AttributeValue::Float(f) => Some(*f),
                        crate::io::AttributeValue::Double(d) => Some(*d as f32),
                        _ => None,
                    })
                    .unwrap_or(-9999.0f32);
                
                let accessor = GenericFieldAccessor {
                    data: data.clone(),
                    longitudes: longitudes.clone(),
                    latitudes: latitudes.clone(),
                    levels: levels_vec.clone(),
                    times: times.clone(),
                    field_name: field_name.to_string(),
                    units,
                    missing_value,
                };
                
                accessors.insert(field_name.to_string(), accessor);
            }
        }
        
        Ok(accessors)
    }
}

impl FieldDataAccess for GenericFieldAccessor {
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
        
        // Check for missing values
        if (value - self.missing_value).abs() < f32::EPSILON {
            return Err(format!("Missing value at indices t={}, k={}, j={}, i={}", time_idx, level_idx, lat_idx, lon_idx));
        }
        
        Ok(value as f64)
    }
    
    fn get_coordinates(&self) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), String> {
        Ok((self.longitudes.clone(), self.latitudes.clone(), self.levels.clone(), self.times.clone()))
    }
}

/// Multi-field accessor for handling multiple meteorological variables
pub struct MultiFieldAccessor {
    pub fields: HashMap<String, GenericFieldAccessor>,
}

impl MultiFieldAccessor {
    /// Create from any DataReader implementation
    pub fn from_data_reader<T: DataReader>(
        reader: &T,
        field_names: &[&str],
    ) -> Result<Self, String> {
        let fields = GenericFieldAccessor::from_data_reader_multi(reader, field_names)?;
        Ok(MultiFieldAccessor { fields })
    }
    
    /// Get accessor for a specific field
    pub fn get_field(&self, field_name: &str) -> Option<&GenericFieldAccessor> {
        self.fields.get(field_name)
    }
    
    /// Get all available field names
    pub fn field_names(&self) -> Vec<&String> {
        self.fields.keys().collect()
    }
    
    /// Interpolate a specific field at a given location and time
    pub fn interpolate_field(
        &self,
        field_name: &str,
        lon: f64,
        lat: f64,
        level: f64,
        time: f64,
    ) -> Result<f64, String> {
        let field_accessor = self.get_field(field_name)
            .ok_or_else(|| format!("Field {} not found", field_name))?;
            
        crate::math::interpolate::interpolate_field_generic(field_accessor, lon, lat, level, time)
    }
}
