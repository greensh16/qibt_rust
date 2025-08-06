use crate::trajectory::TrajectoryPoint;
use crate::data_io::output_trait::{DataWriter, WriteError, TrajectoryMetadata as TraitTrajectoryMetadata};
use chrono::Utc;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::collections::HashMap;

/// NetCDF writer for trajectory output with parallel support
pub struct NetCDFWriter {
    pub file_path: String,
    /// Buffer for rain-cell records during parallel execution
    pub buffer: Arc<Mutex<Vec<BufferedRecord>>>,
    /// Track if file has been created
    pub file_created: bool,
}

/// Buffered trajectory record for parallel writing
#[derive(Debug, Clone)]
pub struct BufferedRecord {
    pub parcel_id: u32,
    pub trajectory: Vec<TrajectoryPoint>,
    pub metadata: TrajectoryMetadata,
}

impl NetCDFWriter {
    /// Create a new NetCDF writer
    pub fn new(file_path: impl AsRef<Path>) -> Self {
        Self {
            file_path: file_path.as_ref().to_string_lossy().to_string(),
            buffer: Arc::new(Mutex::new(Vec::new())),
            file_created: false,
        }
    }

    /// Create output with proper structure and buffered records
    /// Builds file, dimensions, variables, attributes just like Fortran `new_out_file`
    /// During parallel run, buffer each rain-cell record; after join, write sequentially to guarantee correct unlimited dimension order.
    pub fn create_output(
        &mut self,
        _expected_points: usize,
        records: Vec<BufferedRecord>,
    ) -> Result<(), String> {
        // Ensure file is created only once
        if self.file_created {
            return Ok(());
        }

        if records.is_empty() {
            return Err("No trajectory records to write".to_string());
        }

        // Sort records by parcel_id for consistent output order
        let mut sorted_records = records;
        sorted_records.sort_by_key(|r| r.parcel_id);

        println!(
            "Creating NetCDF file: {} with {} trajectories",
            self.file_path,
            sorted_records.len()
        );

        // Find maximum trajectory length for time dimension
        let max_time_steps = sorted_records
            .iter()
            .map(|r| r.trajectory.len())
            .max()
            .unwrap_or(0);

        if max_time_steps == 0 {
            return Err("All trajectories are empty".to_string());
        }

        // Create NetCDF file using the netcdf crate
        let mut file = netcdf::create(&self.file_path)
            .map_err(|e| format!("Failed to create NetCDF file '{}': {}", self.file_path, e))?;

        file.add_unlimited_dimension("time")
            .map_err(|e| format!("Failed to create time dimension: {}", e))?;
        file.add_dimension("trajectory", sorted_records.len())
            .map_err(|e| format!("Failed to create trajectory dimension: {}", e))?;

        // Add global attributes first
        file.add_attribute("title", "Back-trajectory output")
            .map_err(|e| format!("Failed to add title attribute: {}", e))?;
        file.add_attribute("institution", "B-TrIMS Trajectory Model")
            .map_err(|e| format!("Failed to add institution attribute: {}", e))?;
        file.add_attribute("source", "qibt_rust trajectory model")
            .map_err(|e| format!("Failed to add source attribute: {}", e))?;
        file.add_attribute("Conventions", "CF-1.6")
            .map_err(|e| format!("Failed to add Conventions attribute: {}", e))?;

        let creation_time = Utc::now().format("%Y-%m-%dT%H:%M:%S%.3fZ").to_string();
        file.add_attribute("history", format!("Created on {}", creation_time))
            .map_err(|e| format!("Failed to add history attribute: {}", e))?;

        // Create trajectory data variables with proper CF conventions
        let var_specs = [
            ("longitude", "degrees_east", "trajectory longitude"),
            ("latitude", "degrees_north", "trajectory latitude"),
            ("pressure", "Pa", "pressure level"),
            ("temperature", "K", "temperature"),
            ("u_wind", "m/s", "eastward wind component"),
            ("v_wind", "m/s", "northward wind component"),
            ("w_wind", "m/s", "upward wind component"),
        ];

        // Create coordinate variables
        {
            let mut time_var = file
                .add_variable::<f64>("time", &["time"])
                .map_err(|e| format!("Failed to create time variable: {}", e))?;

            // Add attributes to time variable using put_attribute
            time_var
                .put_attribute("units", "hours since trajectory start")
                .map_err(|e| format!("Failed to add time units: {}", e))?;
            time_var
                .put_attribute("long_name", "trajectory time")
                .map_err(|e| format!("Failed to add time long_name: {}", e))?;
        }

        // Create all data variables
        for (name, units, long_name) in &var_specs {
            let mut var = file
                .add_variable::<f64>(name, &["time", "trajectory"])
                .map_err(|e| format!("Failed to create {} variable: {}", name, e))?;

            // Add variable attributes
            var.put_attribute("units", *units)
                .map_err(|e| format!("Failed to add units to {}: {}", name, e))?;
            var.put_attribute("long_name", *long_name)
                .map_err(|e| format!("Failed to add long_name to {}: {}", name, e))?;
            var.put_attribute("_FillValue", -9999.0f64)
                .map_err(|e| format!("Failed to add _FillValue to {}: {}", name, e))?;
        }

        // Write time coordinate variable
        {
            let mut time_var = file
                .variable_mut("time")
                .ok_or_else(|| "Time variable not found".to_string())?;
            let time_values: Vec<f64> = (0..max_time_steps).map(|t| t as f64).collect();
            time_var
                .put_values(&time_values, ..)
                .map_err(|e| format!("Failed to write time values: {}", e))?;
        }

        println!(
            "Writing {} trajectories with up to {} time steps each",
            sorted_records.len(),
            max_time_steps
        );

        // Write trajectory data - organize data sequentially to guarantee correct unlimited dimension order
        for (var_name, _, _) in &var_specs {
            let mut var = file
                .variable_mut(var_name)
                .ok_or_else(|| format!("Variable {} not found", var_name))?;

            // Create data array: time x trajectory
            let mut data = vec![-9999.0f64; max_time_steps * sorted_records.len()];

            for (traj_idx, record) in sorted_records.iter().enumerate() {
                for (time_idx, point) in record.trajectory.iter().enumerate() {
                    let data_idx = time_idx * sorted_records.len() + traj_idx;
                    data[data_idx] = match *var_name {
                        "longitude" => point.longitude,
                        "latitude" => point.latitude,
                        "pressure" => point.pressure,
                        "temperature" => point.temperature,
                        "u_wind" => point.u_wind,
                        "v_wind" => point.v_wind,
                        "w_wind" => point.w_wind,
                        _ => -9999.0,
                    };
                }
            }

            // Write data to NetCDF variable
            var.put_values(&data, (.., ..))
                .map_err(|e| format!("Failed to write {} data: {}", var_name, e))?;
        }

        self.file_created = true;
        println!("Successfully created NetCDF file: {}", self.file_path);
        Ok(())
    }

    /// Buffer each rain-cell record during parallel execution
    pub fn buffer_record(
        &mut self,
        parcel_id: u32,
        trajectory: Vec<TrajectoryPoint>,
        metadata: TrajectoryMetadata,
    ) {
        let record = BufferedRecord {
            parcel_id,
            trajectory,
            metadata,
        };
        let mut buffer = self.buffer.lock().unwrap();
        buffer.push(record);
    }

    /// Write trajectory data to NetCDF file without buffering
    pub fn write_trajectory(&self, trajectory: &[TrajectoryPoint]) -> Result<(), String> {
        // TODO: Implement NetCDF writing using netcdf crate

        println!(
            "Writing {} trajectory points to: {}",
            trajectory.len(),
            self.file_path
        );

        // Placeholder implementation - would write:
        // - Dimensions: time, trajectory_id
        // - Variables: longitude, latitude, pressure, time, temperature, etc.

        // For now, just print summary
        if !trajectory.is_empty() {
            let first = &trajectory[0];
            let last = &trajectory[trajectory.len() - 1];
            println!(
                "  Start: ({:.2}, {:.2}) at {:.2} Pa, time {:.2}",
                first.longitude, first.latitude, first.pressure, first.time
            );
            println!(
                "  End:   ({:.2}, {:.2}) at {:.2} Pa, time {:.2}",
                last.longitude, last.latitude, last.pressure, last.time
            );
        }

        Ok(())
    }

    /// Write multiple trajectories (ensemble or multiple parcels)
    pub fn write_trajectories(&self, trajectories: &[Vec<TrajectoryPoint>]) -> Result<(), String> {
        println!(
            "Writing {} trajectories to: {}",
            trajectories.len(),
            self.file_path
        );

        for (i, traj) in trajectories.iter().enumerate() {
            println!("  Trajectory {}: {} points", i, traj.len());
        }

        // TODO: Implement actual NetCDF writing
        Ok(())
    }

    /// Write trajectory metadata
    pub fn write_metadata(&self, metadata: &TrajectoryMetadata) -> Result<(), String> {
        println!("Writing metadata for trajectory starting at:");
        println!(
            "  Location: ({:.2}, {:.2})",
            metadata.start_longitude, metadata.start_latitude
        );
        println!("  Pressure: {:.2} Pa", metadata.start_pressure);
        println!("  Time: {:.2}", metadata.start_time);

        // TODO: Write as NetCDF global attributes
        Ok(())
    }

    /// Create output file with proper structure
    pub fn create_file(&self, expected_points: usize) -> Result<(), String> {
        println!(
            "Creating NetCDF file: {} (expecting {} points)",
            self.file_path, expected_points
        );

        // TODO: Create NetCDF file with proper dimensions and variables:
        // Dimensions:
        //   - time: unlimited or expected_points
        //   - trajectory_id: number of trajectories
        // Variables:
        //   - time(time): trajectory time coordinates
        //   - longitude(time, trajectory_id): longitude values
        //   - latitude(time, trajectory_id): latitude values
        //   - pressure(time, trajectory_id): pressure values
        //   - temperature(time, trajectory_id): temperature values
        //   - u_wind(time, trajectory_id): u-component wind
        //   - v_wind(time, trajectory_id): v-component wind
        //   - w_wind(time, trajectory_id): w-component wind

        Ok(())
    }
}

/// Metadata for trajectory output
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
}

/// New NetCDF writer implementing the DataWriter trait
pub struct NetCDFTrajectoryWriter {
    file_path: String,
    metadata: Option<TraitTrajectoryMetadata>,
    file_created: bool,
    trajectories: Vec<(u32, Vec<TrajectoryPoint>)>,
}

impl NetCDFTrajectoryWriter {
    pub fn new(file_path: &Path) -> Result<Self, WriteError> {
        Ok(Self {
            file_path: file_path.to_string_lossy().to_string(),
            metadata: None,
            file_created: false,
            trajectories: Vec::new(),
        })
    }
}

impl DataWriter for NetCDFTrajectoryWriter {
    fn create(&mut self, expected_trajectories: usize, expected_time_steps: usize) -> Result<(), WriteError> {
        if self.file_created {
            return Ok(());
        }
        
        println!(
            "Creating NetCDF file: {} (expecting {} trajectories with {} time steps each)",
            self.file_path, expected_trajectories, expected_time_steps
        );
        
        self.file_created = true;
        Ok(())
    }
    
    fn write_trajectory(&mut self, trajectory_id: u32, trajectory: &[TrajectoryPoint]) -> Result<(), WriteError> {
        self.trajectories.push((trajectory_id, trajectory.to_vec()));
        Ok(())
    }
    
    fn write_trajectories(&mut self, trajectories: &[(u32, Vec<TrajectoryPoint>)]) -> Result<(), WriteError> {
        for (id, traj) in trajectories {
            self.trajectories.push((*id, traj.clone()));
        }
        Ok(())
    }
    
    fn set_metadata(&mut self, metadata: &TraitTrajectoryMetadata) -> Result<(), WriteError> {
        self.metadata = Some(metadata.clone());
        Ok(())
    }
    
    fn add_global_attribute(&mut self, _name: &str, _value: &str) -> Result<(), WriteError> {
        // Store in metadata for later writing
        Ok(())
    }
    
    fn close(&mut self) -> Result<(), WriteError> {
        if self.trajectories.is_empty() {
            return Err(WriteError::InvalidData("No trajectories to write".to_string()));
        }
        
        // Actually write the NetCDF file
        let max_time_steps = self.trajectories
            .iter()
            .map(|(_, traj)| traj.len())
            .max()
            .unwrap_or(0);
        
        if max_time_steps == 0 {
            return Err(WriteError::InvalidData("All trajectories are empty".to_string()));
        }
        
        // Create NetCDF file
        let mut file = netcdf::create(&self.file_path)
            .map_err(|e| WriteError::IoError(format!("Failed to create NetCDF file: {}", e)))?;
        
        // Add dimensions
        file.add_unlimited_dimension("time")
            .map_err(|e| WriteError::FormatError(format!("Failed to create time dimension: {}", e)))?;
        file.add_dimension("trajectory", self.trajectories.len())
            .map_err(|e| WriteError::FormatError(format!("Failed to create trajectory dimension: {}", e)))?;
        
        // Add global attributes
        file.add_attribute("title", "Back-trajectory output")
            .map_err(|e| WriteError::FormatError(format!("Failed to add title: {}", e)))?;
        file.add_attribute("institution", "qibt_rust Trajectory Model")
            .map_err(|e| WriteError::FormatError(format!("Failed to add institution: {}", e)))?;
        file.add_attribute("Conventions", "CF-1.6")
            .map_err(|e| WriteError::FormatError(format!("Failed to add Conventions: {}", e)))?;
        
        if let Some(ref metadata) = self.metadata {
            file.add_attribute("creation_time", metadata.creation_time.as_str())
                .map_err(|e| WriteError::FormatError(format!("Failed to add creation_time: {}", e)))?;
            file.add_attribute("meteorological_data_source", metadata.meteorological_data_source.as_str())
                .map_err(|e| WriteError::FormatError(format!("Failed to add data source: {}", e)))?;
        }
        
        // Create variables
        let var_specs = [
            ("longitude", "degrees_east", "trajectory longitude"),
            ("latitude", "degrees_north", "trajectory latitude"),
            ("pressure", "Pa", "pressure level"),
            ("temperature", "K", "temperature"),
            ("u_wind", "m/s", "eastward wind component"),
            ("v_wind", "m/s", "northward wind component"),
            ("w_wind", "m/s", "upward wind component"),
        ];
        
        // Create time variable
        {
            let mut time_var = file
                .add_variable::<f64>("time", &["time"])
                .map_err(|e| WriteError::FormatError(format!("Failed to create time variable: {}", e)))?;
            time_var
                .put_attribute("units", "hours since trajectory start")
                .map_err(|e| WriteError::FormatError(format!("Failed to add time units: {}", e)))?;
            time_var
                .put_attribute("long_name", "trajectory time")
                .map_err(|e| WriteError::FormatError(format!("Failed to add time long_name: {}", e)))?;
        }
        
        // Create data variables
        for (name, units, long_name) in &var_specs {
            let mut var = file
                .add_variable::<f64>(name, &["time", "trajectory"])
                .map_err(|e| WriteError::FormatError(format!("Failed to create {} variable: {}", name, e)))?;
            
            var.put_attribute("units", *units)
                .map_err(|e| WriteError::FormatError(format!("Failed to add units to {}: {}", name, e)))?;
            var.put_attribute("long_name", *long_name)
                .map_err(|e| WriteError::FormatError(format!("Failed to add long_name to {}: {}", name, e)))?;
            var.put_attribute("_FillValue", -9999.0f64)
                .map_err(|e| WriteError::FormatError(format!("Failed to add _FillValue to {}: {}", name, e)))?;
        }
        
        // Write time coordinate
        {
            let mut time_var = file
                .variable_mut("time")
                .ok_or_else(|| WriteError::FormatError("Time variable not found".to_string()))?;
            let time_values: Vec<f64> = (0..max_time_steps).map(|t| t as f64).collect();
            time_var
                .put_values(&time_values, ..)
                .map_err(|e| WriteError::IoError(format!("Failed to write time values: {}", e)))?;
        }
        
        // Write trajectory data
        for (var_name, _, _) in &var_specs {
            let mut var = file
                .variable_mut(var_name)
                .ok_or_else(|| WriteError::FormatError(format!("Variable {} not found", var_name)))?;
            
            let mut data = vec![-9999.0f64; max_time_steps * self.trajectories.len()];
            
            for (traj_idx, (_, trajectory)) in self.trajectories.iter().enumerate() {
                for (time_idx, point) in trajectory.iter().enumerate() {
                    let data_idx = time_idx * self.trajectories.len() + traj_idx;
                    data[data_idx] = match *var_name {
                        "longitude" => point.longitude,
                        "latitude" => point.latitude,
                        "pressure" => point.pressure,
                        "temperature" => point.temperature,
                        "u_wind" => point.u_wind,
                        "v_wind" => point.v_wind,
                        "w_wind" => point.w_wind,
                        _ => -9999.0,
                    };
                }
            }
            
            var.put_values(&data, (.., ..))
                .map_err(|e| WriteError::IoError(format!("Failed to write {} data: {}", var_name, e)))?;
        }
        
        println!("Successfully wrote {} trajectories to NetCDF file: {}", self.trajectories.len(), self.file_path);
        Ok(())
    }
    
    fn get_output_path(&self) -> &str {
        &self.file_path
    }
}

/// Write trajectory to ASCII format (alternative output)
pub fn write_trajectory_ascii(
    file_path: impl AsRef<Path>,
    trajectory: &[TrajectoryPoint],
) -> Result<(), String> {
    use std::fs::File;
    use std::io::Write;

    let mut file = File::create(file_path).map_err(|e| e.to_string())?;

    // Write header
    writeln!(file, "# Back-trajectory output").map_err(|e| e.to_string())?;
    writeln!(file, "# Columns: time(hours) longitude(deg) latitude(deg) pressure(Pa) temperature(K) u_wind(m/s) v_wind(m/s) w_wind(m/s)").map_err(|e| e.to_string())?;

    // Write data
    for point in trajectory {
        writeln!(
            file,
            "{:.6} {:.6} {:.6} {:.2} {:.2} {:.6} {:.6} {:.6}",
            point.time,
            point.longitude,
            point.latitude,
            point.pressure,
            point.temperature,
            point.u_wind,
            point.v_wind,
            point.w_wind
        )
        .map_err(|e| e.to_string())?;
    }

    Ok(())
}
