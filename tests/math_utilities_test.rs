use ndarray::Array1;
use qibt_rust::config::Constants;
use qibt_rust::math::*;

#[test]
fn test_lin_interp() {
    // Test basic linear interpolation
    assert_eq!(lin_interp(1.0, 3.0, 0.5), 2.0);
    assert_eq!(lin_interp(0.0, 10.0, 0.3), 3.0);
    assert_eq!(lin_interp(5.0, 15.0, 0.0), 5.0);
    assert_eq!(lin_interp(5.0, 15.0, 1.0), 15.0);
}

#[test]
fn test_bilin_interp() {
    // Test bilinear interpolation with simple values
    assert_eq!(bilin_interp(1.0, 2.0, 3.0, 4.0, 0.5, 0.5), 2.5);
    assert_eq!(bilin_interp(0.0, 1.0, 2.0, 3.0, 0.0, 0.0), 0.0);
    assert_eq!(bilin_interp(0.0, 1.0, 2.0, 3.0, 1.0, 1.0), 3.0);
}

#[test]
fn test_linear_interpolation() {
    let result = linear_interpolate(0.0, 0.0, 1.0, 10.0, 0.5);
    assert_eq!(result, 5.0);

    let result = linear_interpolate(2.0, 20.0, 4.0, 40.0, 3.0);
    assert_eq!(result, 30.0);
}

#[test]
fn test_bilinear_interpolation() {
    let result = bilinear_interpolate(
        0.0, 0.0, 1.0, 1.0, // Grid bounds
        1.0, 2.0, 3.0, 4.0, // Corner values
        0.5, 0.5, // Interpolation point
    );
    assert_eq!(result, 2.5);
}

#[test]
fn test_lin_interp_array() {
    let v0 = Array1::from(vec![1.0_f64, 2.0_f64, 3.0_f64]);
    let v1 = Array1::from(vec![2.0_f64, 4.0_f64, 6.0_f64]);
    let result = lin_interp_array(&v0, &v1, 0.5_f64);
    let expected = Array1::from(vec![1.5_f64, 3.0_f64, 4.5_f64]);

    for (a, b) in result.iter().zip(expected.iter()) {
        let diff: f64 = *a - *b;
        assert!(diff.abs() < 1e-10);
    }
}

#[test]
fn test_calc_hydrostatic_pressure() {
    let constants = Constants::default();
    let result = calc_hydrostatic_pressure(288.15, 1000.0, &constants);
    // At 1000m height, pressure should be significantly less than surface pressure
    assert!(result < constants.pressure_surface);
    assert!(result > 80000.0); // Should still be reasonable atmospheric pressure
}

#[test]
fn test_calc_pw() {
    let constants = Constants::default();
    let specific_humidity = vec![0.008, 0.006, 0.004]; // kg/kg
    let pressure_levels = vec![100000.0, 90000.0, 80000.0]; // Pa

    let pw = calc_pw(&specific_humidity, &pressure_levels, &constants);

    // Should be positive and reasonable value for precipitable water
    assert!(pw > 0.0);
    assert!(pw < 100.0); // Shouldn't be unreasonably large
}

#[test]
fn test_calc_tpw() {
    let constants = Constants::default();
    let specific_humidity = vec![0.008, 0.006, 0.004];
    let pressure_levels = vec![100000.0, 90000.0, 80000.0];
    let temperature = vec![288.15, 283.15, 278.15]; // K

    let tpw = calc_tpw(
        &specific_humidity,
        &pressure_levels,
        &temperature,
        &constants,
    );

    // Total precipitable water should be positive
    assert!(tpw > 0.0);
    assert!(tpw < 100.0);
}

#[test]
fn test_haversine_distance() {
    // Test distance from equator at 0° longitude to equator at 90° longitude
    let dist = haversine_distance(0.0, 0.0, 0.0, 90.0, 6371000.0);
    let expected = std::f64::consts::PI / 2.0 * 6371000.0; // Quarter of Earth's circumference
    assert!((dist - expected).abs() < 100.0); // Within 100m

    // Test distance between same points (should be 0)
    let dist = haversine_distance(45.0, -100.0, 45.0, -100.0, 6371000.0);
    assert!(dist < 1e-10);
}

#[test]
fn test_pressure_to_height() {
    let constants = Constants::default();
    let height = pressure_to_height(50000.0, 101325.0, &constants);

    // At ~500 hPa, height should be around 5-6 km
    assert!(height > 4000.0);
    assert!(height < 7000.0);
}

#[test]
fn test_height_to_pressure() {
    let constants = Constants::default();
    let pressure = height_to_pressure(5000.0, 101325.0, &constants);

    // At 5 km height, pressure should be around 500-600 hPa
    assert!(pressure > 45000.0);
    assert!(pressure < 65000.0);
}

#[test]
fn test_air_density() {
    let constants = Constants::default();
    let density = air_density(101325.0, 288.15, &constants);

    // Standard atmosphere density at sea level should be ~1.225 kg/m³
    assert!((density - 1.225).abs() < 0.1);
}

#[test]
fn test_potential_temperature() {
    let constants = Constants::default();
    let theta = potential_temperature(288.15, 101325.0, &constants);

    println!(
        "Potential temperature at reference pressure: {}, input: {}",
        theta, 288.15
    );
    // At reference pressure (p0), potential temperature should approximately equal actual temperature
    // But our p0 might be different from the test pressure
    assert!((theta - 288.15).abs() < 5.0); // Allow larger tolerance

    // At lower pressure, potential temperature should be higher
    let theta_high = potential_temperature(260.0, 50000.0, &constants);
    println!("Potential temperature at low pressure: {}", theta_high);
    assert!(theta_high > 260.0);
}

#[test]
fn test_find_grid_indices() {
    let coords = vec![10.0, 20.0, 30.0, 40.0, 50.0];

    // Test interpolation within bounds
    let (i0, i1, weight) = find_grid_indices(&coords, 25.0).unwrap();
    assert_eq!(i0, 1);
    assert_eq!(i1, 2);
    assert!((weight - 0.5).abs() < 1e-10);

    // Test extrapolation below bounds
    let (i0, i1, weight) = find_grid_indices(&coords, 5.0).unwrap();
    assert_eq!(i0, 0);
    assert_eq!(i1, 0);
    assert_eq!(weight, 0.0);

    // Test extrapolation above bounds
    let (i0, i1, weight) = find_grid_indices(&coords, 55.0).unwrap();
    assert_eq!(i0, 4);
    assert_eq!(i1, 4);
    assert_eq!(weight, 0.0);
}

/// Test deterministic arrays for small predictable results
#[test]
fn test_deterministic_interpolation() {
    // Create simple 2x2 test data
    let data = vec![
        vec![vec![vec![1.0, 2.0], vec![3.0, 4.0]]], // time 0, level 0
    ];
    let longitudes = vec![0.0, 1.0];
    let latitudes = vec![0.0, 1.0];
    let levels = vec![1000.0];
    let times = vec![0.0];

    let result = interpolate_meteo_field(
        &data,
        &longitudes,
        &latitudes,
        &levels,
        &times,
        0.5,    // lon
        0.5,    // lat
        1000.0, // level
        0.0,    // time
    )
    .unwrap();

    // Should interpolate to center value: (1+2+3+4)/4 = 2.5
    assert!((result - 2.5).abs() < 1e-10);
}
