pub mod common;
pub mod reader;
pub mod utils;
pub mod writer;

pub use reader::*;
pub use writer::*;
pub mod generic_accessor;

use ndarray::Array4;
use crate::math::interpolate::{FieldDataAccess, ArrayFieldData};

/// Common data structures for meteorological data
#[derive(Debug, Clone)]
pub struct MeteoField {
    /// Field name (e.g., "temperature", "u_wind", "v_wind", "w_wind")
    pub name: String,
    /// 4D data array [time, level, lat, lon] - legacy format
    pub data: Vec<Vec<Vec<Vec<f64>>>>,
    /// Missing value indicator
    pub missing_value: f64,
    /// Units
    pub units: String,
}

/// Modern meteorological field using ndarray
#[derive(Debug, Clone)]
pub struct MeteoFieldArray {
    /// Field name (e.g., "U", "V", "W", "T", "QVAPOR", "RAIN")
    pub name: String,
    /// 4D data array with layout [j, i, k, t] for interpolation-ready access
    pub data: Array4<f32>,
    /// Missing value indicator
    pub missing_value: f32,
    /// Units
    pub units: String,
}

#[derive(Debug, Clone)]
pub struct MeteoGrid {
    /// Longitude coordinates (degrees)
    pub longitudes: Vec<f64>,
    /// Latitude coordinates (degrees)
    pub latitudes: Vec<f64>,
    /// Pressure levels (Pa)
    pub levels: Vec<f64>,
    /// Time coordinates (hours since reference)
    pub times: Vec<f64>,
    /// Reference time for time coordinates
    pub time_reference: String,
}

#[derive(Debug)]
pub struct MeteoData {
    /// Grid definition
    pub grid: MeteoGrid,
    /// Meteorological fields
    pub fields: Vec<MeteoField>,
}

impl Default for MeteoData {
    fn default() -> Self {
        Self::new()
    }
}

impl MeteoData {
    /// Create new empty meteorological data structure
    pub fn new() -> Self {
        Self {
            grid: MeteoGrid {
                longitudes: Vec::new(),
                latitudes: Vec::new(),
                levels: Vec::new(),
                times: Vec::new(),
                time_reference: String::new(),
            },
            fields: Vec::new(),
        }
    }

    /// Get field by name
    pub fn get_field(&self, name: &str) -> Option<&MeteoField> {
        self.fields.iter().find(|field| field.name == name)
    }

    /// Add a new field
    pub fn add_field(&mut self, field: MeteoField) {
        self.fields.push(field);
    }
}
