use crate::config::Constants;
use ndarray::{Array2, Array4};

/// Calculate geopotential height from pressure using hydrostatic equation
pub fn pressure_to_height(pressure: f64, surface_pressure: f64, constants: &Constants) -> f64 {
    let scale_height = constants.r_dry * 288.15 / constants.g; // ~8400m at 288K
    -scale_height * (pressure / surface_pressure).ln()
}

/// Calculate pressure from geopotential height using hydrostatic equation
pub fn height_to_pressure(height: f64, surface_pressure: f64, constants: &Constants) -> f64 {
    let scale_height = constants.r_dry * 288.15 / constants.g;
    surface_pressure * (-height / scale_height).exp()
}

/// Calculate potential temperature from temperature and pressure
pub fn potential_temperature(temperature: f64, pressure: f64, constants: &Constants) -> f64 {
    let kappa = constants.r_dry / constants.cp; // ~0.286
    temperature * (constants.p0 / pressure).powf(kappa)
}

/// Calculate temperature from potential temperature and pressure
pub fn temperature_from_potential(theta: f64, pressure: f64, constants: &Constants) -> f64 {
    let kappa = constants.r_dry / constants.cp;
    theta * (pressure / constants.p0).powf(kappa)
}

/// Calculate air density from pressure and temperature
pub fn air_density(pressure: f64, temperature: f64, constants: &Constants) -> f64 {
    pressure / (constants.r_dry * temperature)
}

/// Calculate hydrostatic pressure
pub fn calc_hydrostatic_pressure(temperature: f64, height: f64, constants: &Constants) -> f64 {
    let scale_height = constants.r_dry * temperature / constants.g;
    constants.pressure_surface * (-height / scale_height).exp()
}

/// Calculate precipitable water vapor (basic version)
pub fn calc_pw(specific_humidity: &[f64], pressure_levels: &[f64], constants: &Constants) -> f64 {
    let mut precipitable_water = 0.0;

    for i in 0..specific_humidity.len() - 1 {
        let dp = pressure_levels[i] - pressure_levels[i + 1]; // pressure difference
        let avg_q = (specific_humidity[i] + specific_humidity[i + 1]) / 2.0; // average specific humidity
        precipitable_water += (avg_q * dp) / constants.g; // integrate using hydrostatic equation
    }

    precipitable_water
}

/// Calculate total precipitable water (comprehensive version)
pub fn calc_tpw(
    specific_humidity: &[f64],
    pressure_levels: &[f64],
    temperature: &[f64],
    constants: &Constants,
) -> f64 {
    let mut total_precipitable_water = 0.0;

    for i in 0..specific_humidity.len() - 1 {
        let dp = pressure_levels[i] - pressure_levels[i + 1];
        let avg_q = (specific_humidity[i] + specific_humidity[i + 1]) / 2.0;
        let avg_temp = (temperature[i] + temperature[i + 1]) / 2.0;

        // Account for temperature dependence
        let correction_factor = 1.0 + (avg_temp - 273.15) * 0.007; // Simple temperature correction
        total_precipitable_water += (avg_q * dp * correction_factor) / constants.g;
    }

    total_precipitable_water
}

/// Simple precipitable water calculation for single layer
pub fn calc_precipitable_water(vapor_pressure: f64, height: f64, constants: &Constants) -> f64 {
    let gas_constant_vapor = constants.rv;
    vapor_pressure * height / (gas_constant_vapor * 273.15) // Assuming standard temp
}

/// Calculate vertical velocity from omega (pressure velocity)
pub fn omega_to_w(omega: f64, pressure: f64, temperature: f64, constants: &Constants) -> f64 {
    // w = -omega * R * T / (p * g)
    -omega * constants.r_dry * temperature / (pressure * constants.g)
}

/// Calculate omega from vertical velocity
pub fn w_to_omega(w: f64, pressure: f64, temperature: f64, constants: &Constants) -> f64 {
    // omega = -w * p * g / (R * T)
    -w * pressure * constants.g / (constants.r_dry * temperature)
}

/// Calculate Coriolis parameter
pub fn coriolis_parameter(latitude: f64) -> f64 {
    let omega_earth = 7.272e-5; // Earth's rotation rate (rad/s)
    2.0 * omega_earth * latitude.to_radians().sin()
}

/// Convert wind components from grid-relative to earth-relative
pub fn grid_to_earth_wind(u_grid: f64, v_grid: f64, rotation_angle: f64) -> (f64, f64) {
    let cos_rot = rotation_angle.cos();
    let sin_rot = rotation_angle.sin();

    let u_earth = u_grid * cos_rot - v_grid * sin_rot;
    let v_earth = u_grid * sin_rot + v_grid * cos_rot;

    (u_earth, v_earth)
}

/// Calculate distance between two geographic points (Haversine formula)
pub fn haversine_distance(lat1: f64, lon1: f64, lat2: f64, lon2: f64, earth_radius: f64) -> f64 {
    let dlat = (lat2 - lat1).to_radians();
    let dlon = (lon2 - lon1).to_radians();
    let lat1_rad = lat1.to_radians();
    let lat2_rad = lat2.to_radians();

    let a =
        (dlat / 2.0).sin().powi(2) + lat1_rad.cos() * lat2_rad.cos() * (dlon / 2.0).sin().powi(2);
    let c = 2.0 * a.sqrt().asin();

    earth_radius * c
}

/// Calculate bearing between two geographic points
pub fn bearing(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    let lat1_rad = lat1.to_radians();
    let lat2_rad = lat2.to_radians();
    let dlon_rad = (lon2 - lon1).to_radians();

    let y = dlon_rad.sin() * lat2_rad.cos();
    let x = lat1_rad.cos() * lat2_rad.sin() - lat1_rad.sin() * lat2_rad.cos() * dlon_rad.cos();

    y.atan2(x).to_degrees()
}

/// Convert latitude/longitude displacement to meters
pub fn deg_to_meters(dlat: f64, dlon: f64, latitude: f64, earth_radius: f64) -> (f64, f64) {
    let lat_rad = latitude.to_radians();
    let dy = dlat.to_radians() * earth_radius;
    let dx = dlon.to_radians() * earth_radius * lat_rad.cos();
    (dx, dy)
}

/// Convert meter displacement to latitude/longitude
pub fn meters_to_deg(dx: f64, dy: f64, latitude: f64, earth_radius: f64) -> (f64, f64) {
    let lat_rad = latitude.to_radians();
    let dlat = dy / earth_radius;
    let dlon = dx / (earth_radius * lat_rad.cos());
    (dlat.to_degrees(), dlon.to_degrees())
}

/// Calculate pressure tendency from vertical motion
pub fn pressure_tendency(w: f64, pressure: f64, temperature: f64, constants: &Constants) -> f64 {
    // dp/dt = -rho * g * w, where rho = p/(R*T)
    let density = air_density(pressure, temperature, constants);
    -density * constants.g * w
}

/// Calculate diabatic heating rate (placeholder implementation)
pub fn diabatic_heating_rate(
    _temperature: f64,
    pressure: f64,
    _humidity: f64,
    _cloud_fraction: f64,
) -> f64 {
    // Simplified heating rate calculation
    // In reality, this would include:
    // - Radiative heating/cooling
    // - Latent heat from phase changes
    // - Sensible heat flux

    // Placeholder: radiative cooling rate
    let cooling_rate = if pressure > 70000.0 {
        -1.0 // K/day in troposphere
    } else {
        -0.5 // K/day in stratosphere
    };

    cooling_rate / 86400.0 // Convert to K/s
}

/// Richardson number for atmospheric stability
pub fn richardson_number(
    du_dz: f64,
    dv_dz: f64,
    dtheta_dz: f64,
    theta: f64,
    constants: &Constants,
) -> f64 {
    let shear_squared = du_dz.powi(2) + dv_dz.powi(2);
    let buoyancy_freq_squared = constants.g / theta * dtheta_dz;

    if shear_squared.abs() < f64::EPSILON {
        f64::INFINITY
    } else {
        buoyancy_freq_squared / shear_squared
    }
}

/// Mixed layer height calculation (simple bulk Richardson method)
pub fn mixed_layer_height(
    surface_theta: f64,
    surface_u: f64,
    surface_v: f64,
    theta_profile: &[f64],
    u_profile: &[f64],
    v_profile: &[f64],
    heights: &[f64],
    critical_ri: f64,
) -> f64 {
    let mut mlt_height = heights[0];

    for i in 1..heights.len() {
        let dtheta = theta_profile[i] - surface_theta;
        let du = u_profile[i] - surface_u;
        let dv = v_profile[i] - surface_v;
        let dz = heights[i] - heights[0];

        let shear_squared = (du.powi(2) + dv.powi(2)) / dz.powi(2);
        let bulk_ri = if shear_squared > f64::EPSILON {
            9.81 * dtheta * dz / (surface_theta * shear_squared * dz.powi(2))
        } else {
            f64::INFINITY
        };

        if bulk_ri > critical_ri {
            mlt_height = heights[i];
            break;
        }
    }

    mlt_height
}

/// Calculate precipitable water from 4D water vapor mixing ratio field (Fortran QIBT style)
/// This matches the calc_pw() function from Fortran QIBT
pub fn calc_pw_4d(
    qvapor: &Array4<f32>,    // Water vapor mixing ratio [j, i, k, t]
    pressure: &Array4<f32>,  // Pressure [j, i, k, t]
    gravity: f32,
) -> Array4<f32> {
    let shape = qvapor.shape();
    let (nj, ni, nk, nt) = (shape[0], shape[1], shape[2], shape[3]);
    let mut pw = Array4::zeros((nj, ni, nk, nt));
    
    // Iterate over horizontal grid points and time
    for t in 0..nt {
        for j in 0..nj {
            for i in 0..ni {
                // Calculate precipitable water at each vertical level
                for k in 1..nk {
                    let dp = pressure[[j, i, k-1, t]] - pressure[[j, i, k, t]];
                    let avg_qv = 0.5 * (qvapor[[j, i, k-1, t]] + qvapor[[j, i, k, t]]);
                    
                    // PW increment = (q_vapor * dp) / g
                    let pw_increment = (avg_qv * dp) / gravity;
                    pw[[j, i, k, t]] = pw[[j, i, k-1, t]] + pw_increment;
                }
            }
        }
    }
    
    pw
}

/// Calculate total precipitable water (TPW) - column integral (Fortran QIBT style)
/// This matches the calc_tpw() function from Fortran QIBT
pub fn calc_tpw_4d(
    qvapor: &Array4<f32>,    // Water vapor mixing ratio [j, i, k, t]
    pressure: &Array4<f32>,  // Pressure [j, i, k, t]
    surface_pressure: &Array4<f32>, // Surface pressure [j, i, 1, t]
    gravity: f32,
) -> Array2<f32> {
    let shape = qvapor.shape();
    let (nj, ni, nk, _nt) = (shape[0], shape[1], shape[2], shape[3]);
    
    // For simplicity, calculate TPW for first time step
    let mut tpw = Array2::zeros((nj, ni));
    
    for j in 0..nj {
        for i in 0..ni {
            let mut column_pw = 0.0;
            
            // Integrate from surface to top
            for k in 1..nk {
                let p_upper = pressure[[j, i, k, 0]];  // Use first time step
                let p_lower = if k == 1 {
                    surface_pressure[[j, i, 0, 0]]
                } else {
                    pressure[[j, i, k-1, 0]]
                };
                
                let dp = p_lower - p_upper;
                let avg_qv = 0.5 * (qvapor[[j, i, k-1.min(nk-1), 0]] + qvapor[[j, i, k, 0]]);
                
                // Add contribution to column integral
                column_pw += (avg_qv * dp) / gravity;
            }
            
            tpw[[j, i]] = column_pw;
        }
    }
    
    tpw
}

/// Calculate actual temperature from perturbation potential temperature (Fortran QIBT style)
/// This matches the calc_actual_temp() function from Fortran QIBT
pub fn calc_actual_temp_4d(
    theta_pert: &Array4<f32>,  // Perturbation potential temperature
    pressure: &Array4<f32>,    // Pressure [Pa]
    base_theta: f32,           // Base potential temperature [K]
    reference_pressure: f32,   // Reference pressure [Pa]
    rd_cp: f32,                // Rd/Cp ratio
) -> Array4<f32> {
    let shape = theta_pert.shape();
    let mut temperature = Array4::zeros((shape[0], shape[1], shape[2], shape[3]));
    
    for ((j, i, k, t), &theta_p) in theta_pert.indexed_iter() {
        let p = pressure[[j, i, k, t]];
        let total_theta = base_theta + theta_p;
        
        // Convert potential temperature to actual temperature
        // T = θ * (P/P0)^(Rd/Cp)
        temperature[[j, i, k, t]] = total_theta * (p / reference_pressure).powf(rd_cp);
    }
    
    temperature
}

/// Calculate quality factor for moisture conservation (Fortran QIBT style)
/// This implements the qfac calculation from Fortran QIBT
pub fn calc_quality_factor(
    parcel_qv: f32,           // Parcel water vapor mixing ratio
    grid_qv: f32,             // Grid point water vapor mixing ratio
    pbl_flag: bool,           // Whether parcel is in PBL
) -> f32 {
    if pbl_flag {
        // In PBL, use full moisture content
        1.0
    } else {
        // Above PBL, apply conservation scaling
        if grid_qv > 0.0 {
            (parcel_qv / grid_qv).min(1.0).max(0.0)
        } else {
            0.0
        }
    }
}

/// Calculate water vapor contribution from a parcel (Fortran QIBT style)
/// This implements the moisture attribution logic from Fortran QIBT
pub fn calc_parcel_moisture_contribution(
    parcel_qv: f32,           // Parcel water vapor mixing ratio [kg/kg]
    grid_tpw: f32,            // Grid total precipitable water [kg/m²]
    quality_factor: f32,      // Quality factor for conservation
    weight: f32,              // Interpolation weight
) -> f32 {
    if grid_tpw > 0.0 {
        (parcel_qv * quality_factor * weight) / grid_tpw
    } else {
        0.0
    }
}

/// Linear interpolation in time between time steps (Fortran QIBT style)
/// This matches the lin_interp_inMM5tsteps() function from Fortran QIBT
pub fn lin_interp_in_timesteps(
    field_t1: &Array4<f32>,   // Field at time t1
    field_t2: &Array4<f32>,   // Field at time t2
    time_weight: f32,         // Weight factor (0 = t1, 1 = t2)
) -> Array4<f32> {
    let shape = field_t1.shape();
    let mut result = Array4::zeros((shape[0], shape[1], shape[2], shape[3]));
    
    for ((j, i, k, t), &val1) in field_t1.indexed_iter() {
        let val2 = field_t2[[j, i, k, t]];
        result[[j, i, k, t]] = val1 * (1.0 - time_weight) + val2 * time_weight;
    }
    
    result
}
