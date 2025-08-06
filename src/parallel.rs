use crate::{
    config::Config,
    data_io::{NetCDFReader, NetCDFWriter},
    trajectory::{LegacyParcel, TrajectoryPoint},
    io::DataReader,
};
use crossbeam_channel::{self, Receiver, Sender};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

/// Parallel trajectory computation using Rayon
pub fn compute_trajectories_parallel(
    config: &Config,
    mut parcels: Vec<LegacyParcel>,
    reader: &NetCDFReader,
) -> Result<Vec<Vec<crate::trajectory::TrajectoryPoint>>, String> {
    println!(
        "Computing {} trajectories in parallel using {} threads",
        parcels.len(),
        rayon::current_num_threads()
    );

    // Process parcels in parallel using Rayon
    let trajectories: Result<Vec<_>, String> = parcels
        .par_iter_mut()
        .enumerate()
        .map(|(i, parcel)| {
            // Each thread gets its own writer (in practice, we'd need thread-safe writers)
            let writer_path = format!("trajectory_{:04}.nc", i);
            let writer = NetCDFWriter::new(&writer_path);

            // Reset parcel to initial conditions
            parcel.reset();

            // Integrate the trajectory
            crate::trajectory::integrate_back_trajectory(config, parcel, reader, &writer)?;

            // Return the trajectory points
            Ok(parcel.trajectory.clone())
        })
        .collect();

    trajectories
}

/// Generic parallel trajectory computation that works with any DataReader
/// This replaces NetCDF-specific assumptions with trait-based access
pub fn compute_trajectories_parallel_generic<T: DataReader + Sync>(
    config: &Config,
    mut parcels: Vec<LegacyParcel>,
    reader: &T,
) -> Result<Vec<Vec<crate::trajectory::TrajectoryPoint>>, String> {
    println!(
        "Computing {} trajectories in parallel using {} threads (generic reader)",
        parcels.len(),
        rayon::current_num_threads()
    );

    // Process parcels in parallel using Rayon
    let trajectories: Result<Vec<_>, String> = parcels
        .par_iter_mut()
        .enumerate()
        .map(|(i, parcel)| {
            // Each thread gets its own writer (in practice, we'd need thread-safe writers)
            let writer_path = format!("trajectory_{:04}.nc", i);
            let writer = NetCDFWriter::new(&writer_path);

            // Reset parcel to initial conditions
            parcel.reset();

            // Integrate the trajectory using generic reader
            crate::trajectory::integrate_back_trajectory_generic(config, parcel, reader, &writer)?;

            // Return the trajectory points
            Ok(parcel.trajectory.clone())
        })
        .collect();

    trajectories
}

/// Generic grid processing for multiple starting locations
pub fn process_grid_parallel_generic<T: DataReader + Sync>(
    config: &Config,
    grid_lons: &[f64],
    grid_lats: &[f64],
    grid_pressures: &[f64],
    reader: &T,
) -> Result<Vec<Vec<crate::trajectory::TrajectoryPoint>>, String> {
    // Create all grid combinations
    let grid_points: Vec<(f64, f64, f64)> = grid_lons
        .iter()
        .flat_map(|&lon| {
            grid_lats.iter().flat_map(move |&lat| {
                grid_pressures
                    .iter()
                    .map(move |&pressure| (lon, lat, pressure))
            })
        })
        .collect();

    println!("Processing {} grid points in parallel (generic reader)", grid_points.len());

    // Process grid points in parallel using Rayon
    let (tx, rx): (Sender<Vec<TrajectoryPoint>>, Receiver<Vec<TrajectoryPoint>>) =
        crossbeam_channel::unbounded();
    let write_mutex = Arc::new(Mutex::new(NetCDFWriter::new("output.nc")));
    let tx_clone = tx.clone();

    // Collect trajectories in a separate thread-safe structure
    let trajectories = Arc::new(Mutex::new(Vec::new()));
    let trajectories_clone = Arc::clone(&trajectories);

    (0..grid_points.len()).into_par_iter().for_each(|idx| {
        let &(lon, lat, pressure) = &grid_points[idx];
        let mut parcel = LegacyParcel::new(idx as u32, lon, lat, pressure, config.start_time);

        // Create a temporary writer for this trajectory
        let writer_path = format!("trajectory_{:06}.nc", idx);
        let writer = NetCDFWriter::new(&writer_path);

        // Integrate trajectory within the thread using generic reader
        if crate::trajectory::integrate_back_trajectory_generic(config, &mut parcel, reader, &writer).is_ok()
        {
            // Store trajectory in the shared collection
            let mut trajectories_guard = trajectories_clone.lock().unwrap();
            trajectories_guard.push(parcel.trajectory.clone());

            // Also send to channel for potential streaming write
            let _ = tx_clone.send(parcel.trajectory.clone());
        }
    });

    // Close the channel sender
    drop(tx);
    drop(tx_clone);

    // Use a consumer thread to write the data from the channel
    std::thread::spawn(move || {
        for trajectory in rx.iter() {
            if let Ok(writer) = write_mutex.lock() {
                let _ = writer.write_trajectory(&trajectory);
            }
        }
    });

    // Return the collected trajectories
    let final_trajectories = trajectories.lock().unwrap().clone();
    Ok(final_trajectories)
}

/// Parallel grid processing for multiple starting locations
pub fn process_grid_parallel(
    config: &Config,
    grid_lons: &[f64],
    grid_lats: &[f64],
    grid_pressures: &[f64],
    reader: &NetCDFReader,
) -> Result<Vec<Vec<crate::trajectory::TrajectoryPoint>>, String> {
    // Create all grid combinations
    let grid_points: Vec<(f64, f64, f64)> = grid_lons
        .iter()
        .flat_map(|&lon| {
            grid_lats.iter().flat_map(move |&lat| {
                grid_pressures
                    .iter()
                    .map(move |&pressure| (lon, lat, pressure))
            })
        })
        .collect();

    println!("Processing {} grid points in parallel", grid_points.len());

    // Process grid points in parallel using Rayon
    let (tx, rx): (Sender<Vec<TrajectoryPoint>>, Receiver<Vec<TrajectoryPoint>>) =
        crossbeam_channel::unbounded();
    let write_mutex = Arc::new(Mutex::new(NetCDFWriter::new("output.nc")));
    let tx_clone = tx.clone();

    // Collect trajectories in a separate thread-safe structure
    let trajectories = Arc::new(Mutex::new(Vec::new()));
    let trajectories_clone = Arc::clone(&trajectories);

    (0..grid_points.len()).into_par_iter().for_each(|idx| {
        let &(lon, lat, pressure) = &grid_points[idx];
        let mut parcel = LegacyParcel::new(idx as u32, lon, lat, pressure, config.start_time);

        // Create a temporary writer for this trajectory
        let writer_path = format!("trajectory_{:06}.nc", idx);
        let writer = NetCDFWriter::new(&writer_path);

        // Integrate trajectory within the thread
        if crate::trajectory::integrate_back_trajectory(config, &mut parcel, reader, &writer).is_ok()
        {
            // Store trajectory in the shared collection
            let mut trajectories_guard = trajectories_clone.lock().unwrap();
            trajectories_guard.push(parcel.trajectory.clone());

            // Also send to channel for potential streaming write
            let _ = tx_clone.send(parcel.trajectory.clone());
        }
    });

    // Close the channel sender
    drop(tx);
    drop(tx_clone);

    // Use a consumer thread to write the data from the channel
    std::thread::spawn(move || {
        for trajectory in rx.iter() {
            if let Ok(writer) = write_mutex.lock() {
                let _ = writer.write_trajectory(&trajectory);
            }
        }
    });

    // Return the collected trajectories
    let final_trajectories = trajectories.lock().unwrap().clone();
    Ok(final_trajectories)
}

/// Parallel processing with custom thread pool
pub fn compute_with_custom_threads(
    config: &Config,
    parcels: Vec<LegacyParcel>,
    reader: &NetCDFReader,
    num_threads: usize,
) -> Result<Vec<Vec<crate::trajectory::TrajectoryPoint>>, String> {
    // Create custom thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|e| format!("Failed to create thread pool: {}", e))?;

    // Use the custom thread pool
    let trajectories = pool.install(|| compute_trajectories_parallel(config, parcels, reader))?;

    Ok(trajectories)
}

/// Chunk-based parallel processing for large datasets
pub fn process_large_dataset(
    config: &Config,
    parcels: Vec<LegacyParcel>,
    reader: &NetCDFReader,
    chunk_size: usize,
) -> Result<Vec<Vec<crate::trajectory::TrajectoryPoint>>, String> {
    let total_parcels = parcels.len();
    let mut all_trajectories = Vec::new();

    // Process in chunks to avoid memory issues
    for (chunk_idx, chunk) in parcels.chunks(chunk_size).enumerate() {
        println!(
            "Processing chunk {} of {} (parcels {}-{})",
            chunk_idx + 1,
            total_parcels.div_ceil(chunk_size),
            chunk_idx * chunk_size + 1,
            (chunk_idx + 1) * chunk_size.min(total_parcels - chunk_idx * chunk_size)
        );

        let chunk_trajectories = compute_trajectories_parallel(config, chunk.to_vec(), reader)?;

        all_trajectories.extend(chunk_trajectories);
    }

    Ok(all_trajectories)
}

/// Parallel quality control for trajectories
pub fn quality_control_parallel(
    trajectories: &[Vec<crate::trajectory::TrajectoryPoint>],
) -> Vec<Vec<String>> {
    trajectories
        .par_iter()
        .map(|trajectory| crate::trajectory::quality_control_trajectory(trajectory))
        .collect()
}

/// Parallel statistics calculation
pub fn calculate_trajectory_statistics_parallel(
    trajectories: &[Vec<crate::trajectory::TrajectoryPoint>],
) -> TrajectoryStatistics {
    use rayon::prelude::*;

    let stats: Vec<_> = trajectories
        .par_iter()
        .filter(|traj| !traj.is_empty())
        .map(|trajectory| {
            let total_distance = calculate_total_distance(trajectory);
            let max_height = trajectory
                .iter()
                .map(|p| {
                    crate::math::physics::pressure_to_height(
                        p.pressure,
                        101325.0,
                        &crate::config::Constants::default(),
                    )
                })
                .fold(f64::NEG_INFINITY, f64::max);

            let min_pressure = trajectory
                .iter()
                .map(|p| p.pressure)
                .fold(f64::INFINITY, f64::min);

            let integration_time = if trajectory.len() > 1 {
                (trajectory.last().unwrap().time - trajectory.first().unwrap().time) * 24.0
            } else {
                0.0
            };

            SingleTrajectoryStats {
                total_distance,
                max_height,
                min_pressure,
                integration_time,
                num_points: trajectory.len(),
            }
        })
        .collect();

    // Aggregate statistics
    let num_trajectories = stats.len();
    if num_trajectories == 0 {
        return TrajectoryStatistics::default();
    }

    let total_distance_sum: f64 = stats.iter().map(|s| s.total_distance).sum();
    let max_height_max: f64 = stats
        .iter()
        .map(|s| s.max_height)
        .fold(f64::NEG_INFINITY, f64::max);
    let min_pressure_min: f64 = stats
        .iter()
        .map(|s| s.min_pressure)
        .fold(f64::INFINITY, f64::min);
    let integration_time_avg: f64 =
        stats.iter().map(|s| s.integration_time).sum::<f64>() / num_trajectories as f64;
    let total_points: usize = stats.iter().map(|s| s.num_points).sum();

    TrajectoryStatistics {
        num_trajectories,
        total_points,
        avg_distance: total_distance_sum / num_trajectories as f64,
        max_height: max_height_max,
        min_pressure: min_pressure_min,
        avg_integration_time: integration_time_avg,
    }
}

/// Calculate total distance for a single trajectory
fn calculate_total_distance(trajectory: &[crate::trajectory::TrajectoryPoint]) -> f64 {
    let mut total_distance = 0.0;

    for i in 1..trajectory.len() {
        let prev = &trajectory[i - 1];
        let curr = &trajectory[i];

        let distance = crate::math::physics::haversine_distance(
            prev.latitude,
            prev.longitude,
            curr.latitude,
            curr.longitude,
            6371000.0, // Earth radius in meters
        );
        total_distance += distance;
    }

    total_distance
}

/// Statistics for a single trajectory
#[derive(Debug, Clone)]
struct SingleTrajectoryStats {
    total_distance: f64,
    max_height: f64,
    min_pressure: f64,
    integration_time: f64,
    num_points: usize,
}

/// Aggregated statistics for all trajectories
#[derive(Debug, Clone)]
pub struct TrajectoryStatistics {
    pub num_trajectories: usize,
    pub total_points: usize,
    pub avg_distance: f64,
    pub max_height: f64,
    pub min_pressure: f64,
    pub avg_integration_time: f64,
}

impl Default for TrajectoryStatistics {
    fn default() -> Self {
        Self {
            num_trajectories: 0,
            total_points: 0,
            avg_distance: 0.0,
            max_height: 0.0,
            min_pressure: 101325.0,
            avg_integration_time: 0.0,
        }
    }
}

impl std::fmt::Display for TrajectoryStatistics {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Trajectory Statistics:")?;
        writeln!(f, "  Number of trajectories: {}", self.num_trajectories)?;
        writeln!(f, "  Total trajectory points: {}", self.total_points)?;
        writeln!(
            f,
            "  Average distance: {:.2} km",
            self.avg_distance / 1000.0
        )?;
        writeln!(f, "  Maximum height: {:.2} m", self.max_height)?;
        writeln!(
            f,
            "  Minimum pressure: {:.2} hPa",
            self.min_pressure / 100.0
        )?;
        writeln!(
            f,
            "  Average integration time: {:.2} hours",
            self.avg_integration_time
        )?;
        Ok(())
    }
}

/// Monitor parallel processing progress
pub struct ProgressMonitor {
    total_tasks: usize,
    completed_tasks: std::sync::atomic::AtomicUsize,
    start_time: std::time::Instant,
}

impl ProgressMonitor {
    pub fn new(total_tasks: usize) -> Self {
        Self {
            total_tasks,
            completed_tasks: std::sync::atomic::AtomicUsize::new(0),
            start_time: std::time::Instant::now(),
        }
    }

    pub fn increment(&self) {
        let completed = self
            .completed_tasks
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed)
            + 1;

        if completed % 100 == 0 || completed == self.total_tasks {
            let elapsed = self.start_time.elapsed();
            let rate = completed as f64 / elapsed.as_secs_f64();
            let eta_seconds = if rate > 0.0 {
                (self.total_tasks - completed) as f64 / rate
            } else {
                0.0
            };

            println!(
                "Progress: {}/{} ({:.1}%) - Rate: {:.1} traj/s - ETA: {:.0}s",
                completed,
                self.total_tasks,
                completed as f64 / self.total_tasks as f64 * 100.0,
                rate,
                eta_seconds
            );
        }
    }
}

/// Parallel trajectory computation with progress monitoring
pub fn compute_trajectories_with_progress(
    config: &Config,
    mut parcels: Vec<LegacyParcel>,
    reader: &NetCDFReader,
) -> Result<Vec<Vec<crate::trajectory::TrajectoryPoint>>, String> {
    let monitor = std::sync::Arc::new(ProgressMonitor::new(parcels.len()));

    let trajectories: Result<Vec<_>, String> = parcels
        .par_iter_mut()
        .enumerate()
        .map(|(i, parcel)| {
            let writer_path = format!("trajectory_{:04}.nc", i);
            let writer = NetCDFWriter::new(&writer_path);

            parcel.reset();
            let result =
                crate::trajectory::integrate_back_trajectory(config, parcel, reader, &writer);

            monitor.increment();

            result.map(|_| parcel.trajectory.clone())
        })
        .collect();

    trajectories
}
