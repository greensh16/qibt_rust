use qibt_rust::config::Config;
use qibt_rust::trajectory::back_traj::implicit_back_traj_w;
use qibt_rust::trajectory::{
    Parcel, parcel_release_height, parcel_release_time,
};

#[test]
fn test_parcel_release_time() {
    let config = Config::default();
    let time = parcel_release_time(&config);
    assert!(time >= config.start_time);
    assert!(time <= config.end_time);
}

#[test]
fn test_parcel_release_height() {
    let config = Config::default();
    let height = parcel_release_height(&config);
    assert!(height >= config.constants.pressure_top);
    assert!(height <= config.constants.pressure_surface);
}

#[test]
fn test_implicit_back_traj_w() {
    let config = Config::default();
    let mut parcel = Parcel {
        lon: 0.0,
        lat: 0.0,
        lev: 0.0,
        pres: 100000.0,
        q: 0.0,
        id: 1,
        time: 2451545.0, // Arbitrary Julian Date
        is_active: true,
    };

    // Wind velocities (u, v, w)
    let u_wind = 10.0;
    let v_wind = 5.0;
    let w_wind = 0.01;

    let dt = 60.0; // 1 minute time step

    // Apply the implicit back trajectory function
    implicit_back_traj_w(&mut parcel, u_wind, v_wind, w_wind, dt as f32, &config);

    // These are arbitrary checks; replace with actual expected values once Fortran benchmarks are available
    assert!(parcel.lon != 0.0);
    assert!(parcel.lat != 0.0);
    assert!(parcel.pres > 100000.0);
}
