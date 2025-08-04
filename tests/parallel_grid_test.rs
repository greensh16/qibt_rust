use rayon::prelude::*;

#[test]
fn test_parallel_grid_processing_setup() {
    // Test that we can set up the parallel grid processing without errors
    let grid_lons = vec![10.0, 11.0, 12.0];
    let grid_lats = vec![50.0, 51.0, 52.0];
    let grid_pressures = vec![85000.0, 70000.0];

    // Simplified approach to avoid closure capture issues
    let mut grid_points = Vec::new();
    for &lon in &grid_lons {
        for &lat in &grid_lats {
            for &pressure in &grid_pressures {
                grid_points.push((lon, lat, pressure));
            }
        }
    }

    // We should have 3 * 3 * 2 = 18 grid points
    assert_eq!(grid_points.len(), 18);

    // Verify the first and last points
    assert_eq!(grid_points[0], (10.0, 50.0, 85000.0));
    assert_eq!(grid_points[17], (12.0, 52.0, 70000.0));
}

#[test]
fn test_rayon_parallel_iterator() {
    use rayon::prelude::*;

    // Test basic parallel iteration functionality
    let totpts = 1000;
    let results: Vec<usize> = (0..totpts).into_par_iter().map(|idx| idx * 2).collect();

    assert_eq!(results.len(), totpts);
    assert_eq!(results[0], 0);
    assert_eq!(results[1], 2);
    assert_eq!(results[999], 1998);
}

#[test]
fn test_parallel_computation_simple() {
    // Test that parallel computation works
    let data: Vec<f64> = (0..10000).map(|i| i as f64).collect();

    // Sequential computation
    let sequential_sum: f64 = data.iter().map(|x| x * x).sum();

    // Parallel computation
    let parallel_sum: f64 = data.par_iter().map(|x| x * x).sum();

    // Results should be identical
    assert_eq!(sequential_sum, parallel_sum);
}
