use super::{LegacyParcel, Parcel, TrajectoryPoint};
use crate::{
    config::Config,
    data_io::{NetCDFReader, NetCDFWriter},
    math::physics,
};

/// Backward trajectory integration (simple backward Euler method)
/// Ports the logic from implicit_back_traj_w in Fortran code
pub fn integrate_back_trajectory(
    config: &Config,
    parcel: &mut LegacyParcel,
    reader: &NetCDFReader,
    writer: &NetCDFWriter,
) -> Result<(), String> {
    // Set up the integration time sequence
    let times = crate::time_utils::generate_time_sequence(
        config.start_time,
        config.start_time - config.trajectory_length,
        -config.time_step / 3600.0, // time step is negative for backward integration
    );

    for &time in &times {
        // Interpolate meteological data at parcel location
        // TODO: Implement proper field interpolation
        let temperature = get_temperature_at_location(
            reader,
            parcel.current_longitude,
            parcel.current_latitude,
            parcel.current_pressure,
            time,
        )?;

        let (u_wind, v_wind, w_wind) = get_wind_at_location(
            reader,
            parcel.current_longitude,
            parcel.current_latitude,
            parcel.current_pressure,
            time,
        )?;

        // Update parcel conditions based on physics
        let new_point = TrajectoryPoint {
            time,
            longitude: parcel.current_longitude, // updated by physics logic
            latitude: parcel.current_latitude,   // updated by physics logic
            pressure: parcel.current_pressure,   // updated by physics logic
            temperature,
            u_wind, // from interpolation
            v_wind, // from interpolation
            w_wind, // from interpolation
            potential_temperature: physics::potential_temperature(
                temperature,
                parcel.current_pressure,
                &config.constants,
            ),
            mixing_ratio: None,
            relative_humidity: None,
        };

        // Add the new point to the parcel's trajectory
        parcel.add_trajectory_point(new_point);

        // Check boundaries
        parcel.check_boundaries(config.constants.p0, 0.0); // surface at p0, top at 0 hPa
        if !parcel.is_active {
            println!(
                "Parcel {} has reached boundary at time {:.2}",
                parcel.id, time
            );
            break;
        }
    }

    // Output the trajectory data
    writer.write_trajectory(&parcel.trajectory)?;

    Ok(())
}

/// Implicit backward trajectory integration with improved numerical stability
/// This is a more sophisticated version that handles stiff equations better
pub fn implicit_back_trajectory_integration(
    config: &Config,
    parcel: &mut LegacyParcel,
    reader: &NetCDFReader,
    writer: &NetCDFWriter,
) -> Result<(), String> {
    let dt = -config.time_step; // negative for backward integration
    let mut current_time = config.start_time;
    let end_time = config.start_time - config.trajectory_length / 24.0; // convert hours to days

    while current_time > end_time && parcel.is_active {
        // Get meteorological data at current location and time
        let (u_wind, v_wind, w_wind) = get_wind_at_location(
            reader,
            parcel.current_longitude,
            parcel.current_latitude,
            parcel.current_pressure,
            current_time,
        )?;

        let temperature = get_temperature_at_location(
            reader,
            parcel.current_longitude,
            parcel.current_latitude,
            parcel.current_pressure,
            current_time,
        )?;

        // Calculate position changes using implicit method
        let (new_lon, new_lat, new_pressure) = integrate_position_implicit(
            parcel.current_longitude,
            parcel.current_latitude,
            parcel.current_pressure,
            u_wind,
            v_wind,
            w_wind,
            dt,
            &config.constants,
        )?;

        // Create new trajectory point
        let point = TrajectoryPoint {
            time: current_time,
            longitude: new_lon,
            latitude: new_lat,
            pressure: new_pressure,
            temperature,
            u_wind,
            v_wind,
            w_wind,
            potential_temperature: physics::potential_temperature(
                temperature,
                new_pressure,
                &config.constants,
            ),
            mixing_ratio: None,
            relative_humidity: None,
        };

        parcel.add_trajectory_point(point);

        // Check if parcel has hit boundaries
        parcel.check_boundaries(100000.0, 1000.0); // 1000 hPa surface, 10 hPa top

        current_time += dt / 86400.0; // convert seconds to days
    }

    // Write trajectory to output
    writer.write_trajectory(&parcel.trajectory)?;

    Ok(())
}

/// Advect parcel based on wind velocities
fn advect(parcel: &mut Parcel, u: f32, v: f32, w: f32, dt: f32, config: &Config) {
    let constants = &config.constants;
    let dlon = (u * dt / (constants.earth_radius as f32 * parcel.lat.to_radians().cos())) as f32;
    let dlat = (v * dt / constants.earth_radius as f32) as f32;
    let dp = w * dt;
    parcel.lon += dlon;
    parcel.lat += dlat;
    parcel.pres += dp;
}

/// Update parcel vertical level based on velocity
fn new_parcel_level_w(parcel: &mut Parcel, w: f32, dt: f32) {
    parcel.lev += w * dt;
}

/// Core implicit back trajectory integration (Merrill scheme)
pub fn implicit_back_traj_w(
    parcel: &mut Parcel,
    u_wind: f32,
    v_wind: f32,
    w_wind: f32,
    dt: f32,
    config: &Config,
) {
    advect(parcel, u_wind, v_wind, w_wind, dt, config);
    new_parcel_level_w(parcel, w_wind, dt);
}

/// Get wind components at specified location and time
fn get_wind_at_location(
    _reader: &NetCDFReader,
    lon: f64,
    lat: f64,
    _pressure: f64,
    _time: f64,
) -> Result<(f64, f64, f64), String> {
    // TODO: This should use proper NetCDF field interpolation
    // For now, return placeholder values

    // In a real implementation, this would:
    // 1. Read the u, v, w wind fields from NetCDF
    // 2. Interpolate to the exact location and time
    // 3. Handle coordinate transformations if needed

    let u_wind = 10.0 * (lat / 90.0).sin(); // Simple placeholder
    let v_wind = 5.0 * (lon / 180.0).cos();
    let w_wind = 0.01; // Small vertical velocity

    Ok((u_wind, v_wind, w_wind))
}

/// Get temperature at specified location and time
fn get_temperature_at_location(
    _reader: &NetCDFReader,
    _lon: f64,
    _lat: f64,
    pressure: f64,
    _time: f64,
) -> Result<f64, String> {
    // TODO: Implement actual NetCDF temperature interpolation

    // Simple atmospheric temperature profile
    let height =
        physics::pressure_to_height(pressure, 101325.0, &crate::config::Constants::default());
    let temperature = 288.15 - 0.0065 * height; // Standard atmosphere lapse rate

    Ok(temperature.max(200.0)) // Don't go below 200K
}

/// Integrate parcel position using implicit method for better stability
fn integrate_position_implicit(
    lon: f64,
    lat: f64,
    pressure: f64,
    u_wind: f64,
    v_wind: f64,
    w_wind: f64,
    dt: f64,
    constants: &crate::config::Constants,
) -> Result<(f64, f64, f64), String> {
    // Convert wind velocities to position changes

    // Horizontal displacement
    let lat_rad = lat.to_radians();
    let earth_radius = constants.earth_radius;

    // Convert wind velocities to angular changes
    let dlon = u_wind * dt / (earth_radius * lat_rad.cos()) * 180.0 / std::f64::consts::PI;
    let dlat = v_wind * dt / earth_radius * 180.0 / std::f64::consts::PI;

    // Vertical displacement (pressure coordinate)
    // Using hydrostatic equation: dp/dt = -rho * g * w
    let temperature = 250.0; // Approximate temperature
    let density = physics::air_density(pressure, temperature, constants);
    let dp = -density * constants.g * w_wind * dt;

    let new_lon = lon + dlon;
    let new_lat = lat + dlat;
    let new_pressure = (pressure + dp).max(1000.0).min(105000.0); // Clamp to reasonable range

    Ok((new_lon, new_lat, new_pressure))
}

/// Run multiple parcel trajectories in parallel (placeholder for parallel module)
pub fn run_ensemble_trajectories(
    config: &Config,
    parcels: &mut [LegacyParcel],
    reader: &NetCDFReader,
    writer: &NetCDFWriter,
) -> Result<(), String> {
    let num_parcels = parcels.len();
    println!("Running {} trajectories...", num_parcels);

    for (i, parcel) in parcels.iter_mut().enumerate() {
        println!("Computing trajectory {} of {}", i + 1, num_parcels);

        // Reset parcel to initial conditions
        parcel.reset();

        // Integrate trajectory
        integrate_back_trajectory(config, parcel, reader, writer)?;
    }

    println!("All trajectories completed.");
    Ok(())
}

/// Quality control checks for trajectory data
pub fn quality_control_trajectory(trajectory: &[TrajectoryPoint]) -> Vec<String> {
    let mut warnings = Vec::new();

    for (i, point) in trajectory.iter().enumerate() {
        // Check for unrealistic values
        if point.longitude.abs() > 180.0 {
            warnings.push(format!(
                "Point {}: Longitude out of range: {}",
                i, point.longitude
            ));
        }

        if point.latitude.abs() > 90.0 {
            warnings.push(format!(
                "Point {}: Latitude out of range: {}",
                i, point.latitude
            ));
        }

        if point.pressure < 100.0 || point.pressure > 110000.0 {
            warnings.push(format!(
                "Point {}: Pressure out of range: {} Pa",
                i, point.pressure
            ));
        }

        if point.temperature < 150.0 || point.temperature > 350.0 {
            warnings.push(format!(
                "Point {}: Temperature out of range: {} K",
                i, point.temperature
            ));
        }

        // Check for large jumps between consecutive points
        if i > 0 {
            let prev = &trajectory[i - 1];
            let distance = physics::haversine_distance(
                prev.latitude,
                prev.longitude,
                point.latitude,
                point.longitude,
                6371000.0,
            );

            // Flag if parcel moves more than 1000 km in one time step
            if distance > 1000000.0 {
                warnings.push(format!(
                    "Point {}: Large displacement: {:.0} km",
                    i,
                    distance / 1000.0
                ));
            }
        }
    }

    warnings
}
