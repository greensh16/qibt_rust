use crate::trajectory::TrajectoryPoint;
use crate::data_io::output_trait::{DataWriter, WriteError, TrajectoryMetadata};
use std::path::Path;
use std::fs::File;
use std::io::Write;

/// ASCII writer for trajectory output
pub struct AsciiTrajectoryWriter {
    file_path: String,
    metadata: Option<TrajectoryMetadata>,
    trajectories: Vec<(u32, Vec<TrajectoryPoint>)>,
}

impl AsciiTrajectoryWriter {
    pub fn new(file_path: &Path) -> Result<Self, WriteError> {
        Ok(Self {
            file_path: file_path.to_string_lossy().to_string(),
            metadata: None,
            trajectories: Vec::new(),
        })
    }
}

impl DataWriter for AsciiTrajectoryWriter {
    fn create(&mut self, expected_trajectories: usize, expected_time_steps: usize) -> Result<(), WriteError> {
        println!(
            "Creating ASCII file: {} (expecting {} trajectories with {} time steps each)",
            self.file_path, expected_trajectories, expected_time_steps
        );
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
    
    fn set_metadata(&mut self, metadata: &TrajectoryMetadata) -> Result<(), WriteError> {
        self.metadata = Some(metadata.clone());
        Ok(())
    }
    
    fn add_global_attribute(&mut self, _name: &str, _value: &str) -> Result<(), WriteError> {
        // ASCII files don't support structured attributes
        Ok(())
    }
    
    fn close(&mut self) -> Result<(), WriteError> {
        if self.trajectories.is_empty() {
            return Err(WriteError::InvalidData("No trajectories to write".to_string()));
        }
        
        let mut file = File::create(&self.file_path)
            .map_err(|e| WriteError::IoError(format!("Failed to create ASCII file: {}", e)))?;
        
        // Write header
        writeln!(file, "# Back-trajectory output")
            .map_err(|e| WriteError::IoError(e.to_string()))?;
        
        if let Some(ref metadata) = self.metadata {
            writeln!(file, "# Creation time: {}", metadata.creation_time)
                .map_err(|e| WriteError::IoError(e.to_string()))?;
            writeln!(file, "# Data source: {}", metadata.meteorological_data_source)
                .map_err(|e| WriteError::IoError(e.to_string()))?;
            writeln!(file, "# Start location: ({:.6}, {:.6})", metadata.start_longitude, metadata.start_latitude)
                .map_err(|e| WriteError::IoError(e.to_string()))?;
            writeln!(file, "# Start pressure: {:.2} Pa", metadata.start_pressure)
                .map_err(|e| WriteError::IoError(e.to_string()))?;
        }
        
        writeln!(file, "# Columns: trajectory_id time(hours) longitude(deg) latitude(deg) pressure(Pa) temperature(K) u_wind(m/s) v_wind(m/s) w_wind(m/s)")
            .map_err(|e| WriteError::IoError(e.to_string()))?;
        
        // Write trajectory data
        for (traj_id, trajectory) in &self.trajectories {
            for point in trajectory {
                writeln!(
                    file,
                    "{} {:.6} {:.6} {:.6} {:.2} {:.2} {:.6} {:.6} {:.6}",
                    traj_id,
                    point.time,
                    point.longitude,
                    point.latitude,
                    point.pressure,
                    point.temperature,
                    point.u_wind,
                    point.v_wind,
                    point.w_wind
                )
                .map_err(|e| WriteError::IoError(e.to_string()))?;
            }
        }
        
        println!("Successfully wrote {} trajectories to ASCII file: {}", self.trajectories.len(), self.file_path);
        Ok(())
    }
    
    fn get_output_path(&self) -> &str {
        &self.file_path
    }
}
