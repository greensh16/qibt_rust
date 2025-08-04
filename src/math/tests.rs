use super::interpolate::*;
use super::physics::*;
use crate::config::Constants;

#[test]
fn test_lin_interp() {
    assert_eq!(lin_interp(1.0, 3.0, 0.5), 2.0);
}

#[test]
fn test_bilin_interp() {
    assert_eq!(bilin_interp(1.0, 2.0, 3.0, 4.0, 0.5, 0.5), 2.5);
}

#[test]
fn test_calc_hydrostatic_pressure() {
    let constants = Constants::default();
    let result = calc_hydrostatic_pressure(288.15, 1000.0, &constants);
    println!(
        "Hydrostatic pressure result: {}, expected around: {}",
        result, 101325.0
    );
    // At 1000m height with standard temp, pressure should be ~90 kPa, not 101 kPa
    assert!(result < constants.pressure_surface);
    assert!(result > 80000.0);
}

#[test]
fn test_calc_pw() {
    let constants = Constants::default();
    let specific_humidity = vec![0.007, 0.006];
    let pressure_levels = vec![100000.0, 90000.0]; // Convert to Pa from hPa
    let pw = calc_pw(&specific_humidity, &pressure_levels, &constants);
    println!("Precipitable water result: {}", pw);
    assert!(pw > 0.0);
    assert!(pw < 100.0); // Reasonable upper bound
}

#[test]
fn test_haversine_distance() {
    let dist = haversine_distance(0.0, 0.0, 0.0, 90.0, 6371000.0);
    assert!((dist - 10007557.0).abs() < 1000.0);
}

#[test]
fn test_linear_interpolation() {
    let result = linear_interpolate(0.0, 0.0, 1.0, 10.0, 0.5);
    assert_eq!(result, 5.0);
}
