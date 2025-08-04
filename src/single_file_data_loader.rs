use chrono::{DateTime, Utc, TimeZone};
use std::sync::{Arc, Mutex};
use std::collections::HashMap;
use netcdf::File;
use ndarray::Array4;
use crate::time_series_netcdf_reader::TimeSeriesNetCDFReader;

struct SingleFileDataLoader {
    reader: TimeSeriesNetCDFReader,
    time: Vec<DateTime<Utc>>,
}

impl SingleFileDataLoader {
    fn new(file_path: &str) -> Result<Self, String> {
        let file = File::open(file_path).map_err(|e| e.to_string())?;
        let time_var = file.variable("time").ok_or("Missing time variable")?;
        let time_data: Vec<f64> = time_var.get_values::<f64>(..).map_err(|e| e.to_string())?;
        let time: Vec<DateTime<Utc>> = time_data.into_iter().map(|julian| {
            // Convert Julian date to DateTime<Utc>
            Utc.timestamp_opt((julian - 2440587.5) * 86400.0, 0).unwrap()
        }).collect();

        let reader = TimeSeriesNetCDFReader::new(file);

        Ok(SingleFileDataLoader {
            reader,
            time,
        })
    }

    fn get_variable_slice(&self, var: &str, datetime: DateTime<Utc>) -> Result<Array4<f32>, String> {
        // Find the nearest time index
        let time_index = self.time.binary_search(&datetime).unwrap_or_else(|x| x); // nearest index
        
        // Use TimeSeriesNetCDFReader to get the variable array for a single time step
        let time_range = (time_index, time_index + 1);
        
        self.reader.read_variable_array(time_range, var)
            .map_err(|e| format!("Reader error: {}", e))
    }
}

fn main() {
    // Example usage
    let loader = SingleFileDataLoader::new("/path/to/netcdf").expect("Failed to open");
    let time = Utc::now(); // Example time
    
    match loader.get_variable_slice("temperature", time) {
        Ok(slice) => println!("Slice shape: {:?}", slice.shape()),
        Err(e) => eprintln!("Error: {}", e),
    }
}
