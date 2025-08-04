use crate::config::Config;
use crate::data_io::{NetCDFReader, NetCDFWriter};
use crate::parallel::{
    TrajectoryStatistics, calculate_trajectory_statistics_parallel, compute_with_custom_threads,
};
use crate::trajectory::LegacyParcel;
use std::fs;
use std::path::Path;
use std::time::{Duration, Instant};

/// Benchmark results for a specific thread count
#[derive(Debug, Clone)]
pub struct BenchmarkResult {
    pub thread_count: usize,
    pub duration: Duration,
    pub parcels_per_second: f64,
    pub memory_usage_mb: f64,
    pub statistics: TrajectoryStatistics,
}

/// Complete benchmark suite results
#[derive(Debug, Clone)]
pub struct BenchmarkSuite {
    pub results: Vec<BenchmarkResult>,
    pub input_file: String,
    pub parcel_count: usize,
    pub trajectory_length: f64,
}

impl BenchmarkSuite {
    /// Run comprehensive benchmark suite
    pub fn run_suite(
        input_file: &str,
        output_dir: &str,
        parcel_count: usize,
        thread_counts: &[usize],
        enable_flamegraph: bool,
    ) -> Result<Self, String> {
        // Create output directory
        fs::create_dir_all(output_dir)
            .map_err(|e| format!("Failed to create output directory: {}", e))?;

        // Check input file exists
        if !Path::new(input_file).exists() {
            return Err(format!("Input file does not exist: {}", input_file));
        }

        println!("Starting benchmark suite:");
        println!("  Input file: {}", input_file);
        println!("  Parcel count: {}", parcel_count);
        println!("  Thread counts: {:?}", thread_counts);
        println!("  Output directory: {}", output_dir);
        println!("  Flamegraph enabled: {}", enable_flamegraph);

        let mut results = Vec::new();

        for &thread_count in thread_counts {
            println!("\\nBenchmarking with {} threads...", thread_count);

            let result = Self::benchmark_single_configuration(
                input_file,
                output_dir,
                parcel_count,
                thread_count,
                enable_flamegraph,
            )?;

            println!("  Duration: {:.2}s", result.duration.as_secs_f64());
            println!("  Parcels/sec: {:.2}", result.parcels_per_second);
            println!("  Memory usage: {:.2} MB", result.memory_usage_mb);

            results.push(result);
        }

        let suite = BenchmarkSuite {
            results,
            input_file: input_file.to_string(),
            parcel_count,
            trajectory_length: 48.0, // Default 2 days
        };

        // Write benchmark report
        suite.write_report(output_dir)?;

        Ok(suite)
    }

    /// Benchmark a single thread configuration
    fn benchmark_single_configuration(
        input_file: &str,
        output_dir: &str,
        parcel_count: usize,
        thread_count: usize,
        _enable_flamegraph: bool,
    ) -> Result<BenchmarkResult, String> {
        // Create configuration
        let config = Config {
            input_path: Path::new(input_file).to_path_buf(),
            num_parcels: parcel_count,
            num_threads: thread_count,
            trajectory_length: 48.0,
            time_step: 600.0,
            ..Default::default()
        };

        // Create initial parcels
        let parcels: Vec<LegacyParcel> = (0..parcel_count)
            .map(|i| {
                let lat = 40.0 + (i as f64 % 10.0);
                let lon = -100.0 + (i as f64 % 10.0);
                LegacyParcel::new(i as u32, lon, lat, 50000.0, config.start_time)
            })
            .collect();

        // Create mock reader (placeholder for actual NetCDF reader)
        let reader = NetCDFReader::new(input_file);

        // Measure memory before computation
        let memory_before = get_memory_usage_mb();

        // Perform the benchmark
        let start_time = Instant::now();

        let trajectories = compute_with_custom_threads(&config, parcels, &reader, thread_count)
            .map_err(|e| format!("Benchmark computation failed: {}", e))?;

        let duration = start_time.elapsed();

        // Measure memory after computation
        let memory_after = get_memory_usage_mb();
        let memory_usage_mb = memory_after - memory_before;

        // Calculate statistics
        let statistics = calculate_trajectory_statistics_parallel(&trajectories);

        // Calculate performance metrics
        let parcels_per_second = parcel_count as f64 / duration.as_secs_f64();

        // Write output for this benchmark
        let output_file = format!("{}/benchmark_{}threads.nc", output_dir, thread_count);
        let writer = NetCDFWriter::new(&output_file);
        for trajectory in &trajectories {
            let _ = writer.write_trajectory(trajectory);
        }

        Ok(BenchmarkResult {
            thread_count,
            duration,
            parcels_per_second,
            memory_usage_mb,
            statistics,
        })
    }

    /// Write comprehensive benchmark report
    pub fn write_report(&self, output_dir: &str) -> Result<(), String> {
        let report_path = format!("{}/benchmark_report.txt", output_dir);
        let mut report = String::new();

        report.push_str("=== QIBT Rust Benchmark Report ===\\n\\n");
        report.push_str(&format!("Input file: {}\\n", self.input_file));
        report.push_str(&format!("Parcel count: {}\\n", self.parcel_count));
        report.push_str(&format!(
            "Trajectory length: {:.1} hours\\n",
            self.trajectory_length
        ));
        report.push_str("\\n");

        // Performance summary table
        report.push_str("Performance Summary:\\n");
        report.push_str("Threads | Duration (s) | Parcels/sec | Memory (MB) | Speedup\\n");
        report.push_str("--------|-------------|-------------|-------------|--------\\n");

        let baseline_duration = self
            .results
            .first()
            .map(|r| r.duration.as_secs_f64())
            .unwrap_or(1.0);

        for result in &self.results {
            let speedup = baseline_duration / result.duration.as_secs_f64();
            report.push_str(&format!(
                "{:7} | {:11.2} | {:11.2} | {:11.2} | {:6.2}x\\n",
                result.thread_count,
                result.duration.as_secs_f64(),
                result.parcels_per_second,
                result.memory_usage_mb,
                speedup
            ));
        }

        report.push_str("\\n");

        // Detailed statistics for each configuration
        for result in &self.results {
            report.push_str(&format!(
                "\\n=== {} Threads Detailed Stats ===\\n",
                result.thread_count
            ));
            report.push_str(&format!("{}", result.statistics));
        }

        // Performance analysis
        report.push_str("\\n=== Performance Analysis ===\\n");
        if self.results.len() >= 2 {
            let single_thread = &self.results[0];
            let best_result = self
                .results
                .iter()
                .min_by(|a, b| a.duration.partial_cmp(&b.duration).unwrap())
                .unwrap();

            let max_speedup =
                single_thread.duration.as_secs_f64() / best_result.duration.as_secs_f64();
            report.push_str(&format!(
                "Maximum speedup: {:.2}x with {} threads\\n",
                max_speedup, best_result.thread_count
            ));

            // Scaling efficiency
            let efficiency = max_speedup / best_result.thread_count as f64 * 100.0;
            report.push_str(&format!("Parallel efficiency: {:.1}%\\n", efficiency));
        }

        // Recommendations
        report.push_str("\\n=== Recommendations ===\\n");
        if let Some(optimal) = self.find_optimal_thread_count() {
            report.push_str(&format!("Optimal thread count: {}\\n", optimal));
        }
        report.push_str("Consider using f32 precision for additional performance gains.\\n");

        std::fs::write(&report_path, report)
            .map_err(|e| format!("Failed to write benchmark report: {}", e))?;

        println!("\\nBenchmark report written to: {}", report_path);

        Ok(())
    }

    /// Find the optimal thread count based on performance
    pub fn find_optimal_thread_count(&self) -> Option<usize> {
        self.results
            .iter()
            .min_by(|a, b| a.duration.partial_cmp(&b.duration).unwrap())
            .map(|result| result.thread_count)
    }

    /// Get scaling efficiency for a given thread count
    pub fn get_scaling_efficiency(&self, thread_count: usize) -> Option<f64> {
        let single_thread_time = self.results.first()?.duration.as_secs_f64();

        let target_result = self
            .results
            .iter()
            .find(|r| r.thread_count == thread_count)?;

        let speedup = single_thread_time / target_result.duration.as_secs_f64();
        Some(speedup / thread_count as f64 * 100.0)
    }
}

/// Get current memory usage in MB (placeholder implementation)
fn get_memory_usage_mb() -> f64 {
    // This is a placeholder - in a real implementation, you'd use system-specific
    // memory monitoring or a crate like `sys-info`
    0.0
}

/// Run flamegraph profiling for a specific configuration
pub fn run_flamegraph_profiling(
    input_file: &str,
    output_dir: &str,
    parcel_count: usize,
    thread_count: usize,
) -> Result<(), String> {
    println!(
        "Running flamegraph profiling with {} threads...",
        thread_count
    );

    // This would integrate with cargo flamegraph
    // For now, just run a normal benchmark and suggest manual flamegraph usage
    println!("To run flamegraph profiling manually:");
    println!("  cargo install flamegraph");
    println!("  sudo cargo flamegraph --bin qibt_rust -- benchmark \\");
    println!("    --input {} \\", input_file);
    println!("    --output-dir {} \\", output_dir);
    println!("    --parcels {} \\", parcel_count);
    println!("    --thread-counts {}", thread_count);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_benchmark_result_creation() {
        let result = BenchmarkResult {
            thread_count: 8,
            duration: Duration::from_secs(10),
            parcels_per_second: 100.0,
            memory_usage_mb: 50.0,
            statistics: TrajectoryStatistics::default(),
        };

        assert_eq!(result.thread_count, 8);
        assert_eq!(result.duration, Duration::from_secs(10));
        assert_eq!(result.parcels_per_second, 100.0);
    }
}
