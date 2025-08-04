use crate::config::Config;
use rand::Rng;

/// Air parcel for trajectory physics (matching Merrill scheme)
/// Using f32 to match Fortran precision requirements
#[derive(Debug, Clone, PartialEq)]
pub struct Parcel {
    /// Longitude in degrees
    pub lon: f32,
    /// Latitude in degrees
    pub lat: f32,
    /// Vertical level (model level or pressure)
    pub lev: f32,
    /// Pressure in Pa
    pub pres: f32,
    /// Specific humidity (mixing ratio) in kg/kg
    pub q: f32,
    /// Parcel ID for tracking
    pub id: u32,
    /// Time (Julian day)
    pub time: f64,
    /// Status flags
    pub is_active: bool,
}

/// Represents a single point along a trajectory (original structure preserved for compatibility)
#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryPoint {
    /// Time (Julian day)
    pub time: f64,
    /// Longitude (degrees)
    pub longitude: f64,
    /// Latitude (degrees)
    pub latitude: f64,
    /// Pressure (Pa)
    pub pressure: f64,
    /// Temperature (K)
    pub temperature: f64,
    /// U-component wind velocity (m/s)
    pub u_wind: f64,
    /// V-component wind velocity (m/s)
    pub v_wind: f64,
    /// Vertical wind velocity (m/s)
    pub w_wind: f64,
    /// Potential temperature (K)
    pub potential_temperature: f64,
    /// Mixing ratio (kg/kg) - optional
    pub mixing_ratio: Option<f64>,
    /// Relative humidity (%) - optional
    pub relative_humidity: Option<f64>,
}

impl TrajectoryPoint {
    /// Create a new trajectory point
    pub fn new(time: f64, longitude: f64, latitude: f64, pressure: f64) -> Self {
        Self {
            time,
            longitude,
            latitude,
            pressure,
            temperature: 0.0,
            u_wind: 0.0,
            v_wind: 0.0,
            w_wind: 0.0,
            potential_temperature: 0.0,
            mixing_ratio: None,
            relative_humidity: None,
        }
    }

    /// Get position as (longitude, latitude, pressure) tuple
    pub fn position(&self) -> (f64, f64, f64) {
        (self.longitude, self.latitude, self.pressure)
    }

    /// Get wind vector as (u, v, w) tuple
    pub fn wind_vector(&self) -> (f64, f64, f64) {
        (self.u_wind, self.v_wind, self.w_wind)
    }

    /// Calculate wind speed (horizontal)
    pub fn wind_speed(&self) -> f64 {
        (self.u_wind.powi(2) + self.v_wind.powi(2)).sqrt()
    }

    /// Calculate wind direction (degrees from north, meteorological convention)
    pub fn wind_direction(&self) -> f64 {
        let dir = self.v_wind.atan2(self.u_wind).to_degrees();
        // Convert to meteorological convention (direction FROM which wind blows)
        (270.0 - dir).rem_euclid(360.0)
    }
}

/// Legacy air parcel with trajectory tracking capabilities (for compatibility)
#[derive(Debug, Clone)]
pub struct LegacyParcel {
    /// Unique identifier for this parcel
    pub id: u32,
    /// Starting position
    pub start_longitude: f64,
    pub start_latitude: f64,
    pub start_pressure: f64,
    pub start_time: f64,
    /// Current position
    pub current_longitude: f64,
    pub current_latitude: f64,
    pub current_pressure: f64,
    pub current_time: f64,
    /// Trajectory history
    pub trajectory: Vec<TrajectoryPoint>,
    /// Parcel properties
    pub mass: f64, // kg
    pub density: f64, // kg/m³
    pub volume: f64,  // m³
    /// Status flags
    pub is_active: bool,
    pub has_reached_surface: bool,
    pub has_reached_top: bool,
}

impl Parcel {
    /// Create a new air parcel at specified location and time
    pub fn new(id: u32, longitude: f64, latitude: f64, pressure: f64, time: f64) -> Self {
        Self {
            id,
            lon: longitude as f32,
            lat: latitude as f32,
            lev: 0.0, // Initialize to 0
            pres: pressure as f32,
            q: 0.001, // Initialize to 1 g/kg
            time,
            is_active: true,
        }
    }

    /// Update parcel position and properties
    pub fn update(&mut self, longitude: f64, latitude: f64, pressure: f64, time: f64) {
        self.lon = longitude as f32;
        self.lat = latitude as f32;
        self.pres = pressure as f32;
        self.time = time;
    }

    /// Get position as (longitude, latitude, pressure) tuple
    pub fn position(&self) -> (f64, f64, f64) {
        (self.lon as f64, self.lat as f64, self.pres as f64)
    }

    /// Check if parcel has left the model domain
    pub fn is_out_of_bounds(&self, lon_min: f64, lon_max: f64, lat_min: f64, lat_max: f64) -> bool {
        (self.lon as f64) < lon_min
            || (self.lon as f64) > lon_max
            || (self.lat as f64) < lat_min
            || (self.lat as f64) > lat_max
    }

    /// Check if parcel has reached surface or top boundary
    pub fn check_boundaries(&mut self, surface_pressure: f64, top_pressure: f64) {
        if (self.pres as f64) >= surface_pressure {
            self.is_active = false;
        }

        if (self.pres as f64) <= top_pressure {
            self.is_active = false;
        }
    }
}

impl LegacyParcel {
    /// Create a new air parcel at specified location and time
    pub fn new(id: u32, longitude: f64, latitude: f64, pressure: f64, time: f64) -> Self {
        Self {
            id,
            start_longitude: longitude,
            start_latitude: latitude,
            start_pressure: pressure,
            start_time: time,
            current_longitude: longitude,
            current_latitude: latitude,
            current_pressure: pressure,
            current_time: time,
            trajectory: Vec::new(),
            mass: 1.0,
            density: 1.225, // Standard air density at sea level
            volume: 1.0,
            is_active: true,
            has_reached_surface: false,
            has_reached_top: false,
        }
    }

    /// Add a new point to the trajectory
    pub fn add_trajectory_point(&mut self, point: TrajectoryPoint) {
        self.current_longitude = point.longitude;
        self.current_latitude = point.latitude;
        self.current_pressure = point.pressure;
        self.current_time = point.time;
        self.trajectory.push(point);
    }

    /// Update parcel position
    pub fn update_position(&mut self, longitude: f64, latitude: f64, pressure: f64, time: f64) {
        self.current_longitude = longitude;
        self.current_latitude = latitude;
        self.current_pressure = pressure;
        self.current_time = time;
    }

    /// Get the most recent trajectory point
    pub fn last_point(&self) -> Option<&TrajectoryPoint> {
        self.trajectory.last()
    }

    /// Calculate total distance traveled
    pub fn total_distance(&self) -> f64 {
        use crate::config::Constants;
        use crate::math::physics::haversine_distance;

        let constants = Constants::default();
        let mut total_distance = 0.0;

        for i in 1..self.trajectory.len() {
            let prev = &self.trajectory[i - 1];
            let curr = &self.trajectory[i];

            let distance = haversine_distance(
                prev.latitude,
                prev.longitude,
                curr.latitude,
                curr.longitude,
                constants.earth_radius,
            );
            total_distance += distance;
        }

        total_distance
    }

    /// Calculate integration time (hours)
    pub fn integration_time(&self) -> f64 {
        if self.trajectory.is_empty() {
            0.0
        } else {
            (self.current_time - self.start_time) * 24.0
        }
    }

    /// Check if parcel has left the model domain
    pub fn is_out_of_bounds(&self, lon_min: f64, lon_max: f64, lat_min: f64, lat_max: f64) -> bool {
        self.current_longitude < lon_min
            || self.current_longitude > lon_max
            || self.current_latitude < lat_min
            || self.current_latitude > lat_max
    }

    /// Check if parcel has reached surface or top boundary
    pub fn check_boundaries(&mut self, surface_pressure: f64, top_pressure: f64) {
        if self.current_pressure >= surface_pressure {
            self.has_reached_surface = true;
            self.is_active = false;
        }

        if self.current_pressure <= top_pressure {
            self.has_reached_top = true;
            self.is_active = false;
        }
    }

    /// Reset parcel to initial conditions
    pub fn reset(&mut self) {
        self.current_longitude = self.start_longitude;
        self.current_latitude = self.start_latitude;
        self.current_pressure = self.start_pressure;
        self.current_time = self.start_time;
        self.trajectory.clear();
        self.is_active = true;
        self.has_reached_surface = false;
        self.has_reached_top = false;
    }
}

/// Parcel release configuration
#[derive(Debug, Clone)]
pub struct ReleaseConfig {
    /// Release location bounds
    pub lon_min: f64,
    pub lon_max: f64,
    pub lat_min: f64,
    pub lat_max: f64,
    /// Pressure level bounds
    pub pressure_min: f64,
    pub pressure_max: f64,
    /// Time bounds
    pub time_start: f64,
    pub time_end: f64,
    /// Grid spacing
    pub lon_spacing: f64,
    pub lat_spacing: f64,
    pub pressure_spacing: f64,
    pub time_spacing: f64,
    /// Release pattern
    pub release_pattern: ReleasePattern,
}

/// Different patterns for releasing parcels
#[derive(Debug, Clone)]
pub enum ReleasePattern {
    /// Regular grid in space and time
    Grid,
    /// Random distribution
    Random { count: usize },
    /// Single point release
    Point,
    /// Line release (e.g., along flight track)
    Line { waypoints: Vec<(f64, f64, f64)> },
    /// Area source (uniform over specified area)
    Area,
}

/// Generate parcels according to release configuration
pub fn generate_parcels(config: &ReleaseConfig) -> Vec<Parcel> {
    match &config.release_pattern {
        ReleasePattern::Grid => generate_grid_parcels(config),
        ReleasePattern::Random { count } => generate_random_parcels(config, *count),
        ReleasePattern::Point => generate_point_parcels(config),
        ReleasePattern::Line { waypoints } => generate_line_parcels(config, waypoints),
        ReleasePattern::Area => generate_area_parcels(config),
    }
}

/// Generate parcels on a regular grid
fn generate_grid_parcels(config: &ReleaseConfig) -> Vec<Parcel> {
    let mut parcels = Vec::new();
    let mut id = 0;

    let lon_steps = ((config.lon_max - config.lon_min) / config.lon_spacing).ceil() as usize + 1;
    let lat_steps = ((config.lat_max - config.lat_min) / config.lat_spacing).ceil() as usize + 1;
    let pressure_steps =
        ((config.pressure_max - config.pressure_min) / config.pressure_spacing).ceil() as usize + 1;
    let time_steps =
        ((config.time_end - config.time_start) / config.time_spacing).ceil() as usize + 1;

    for it in 0..time_steps {
        let time = config.time_start + it as f64 * config.time_spacing;

        for ip in 0..pressure_steps {
            let pressure = config.pressure_min + ip as f64 * config.pressure_spacing;

            for ilat in 0..lat_steps {
                let lat = config.lat_min + ilat as f64 * config.lat_spacing;

                for ilon in 0..lon_steps {
                    let lon = config.lon_min + ilon as f64 * config.lon_spacing;

                    parcels.push(Parcel::new(id, lon, lat, pressure, time));
                    id += 1;
                }
            }
        }
    }

    parcels
}

/// Generate randomly distributed parcels
/// Random precipitation-weighted sampling for release time
pub fn parcel_release_time(config: &Config) -> f64 {
    let mut rng = rand::rng();
    let precip_weight = rng.random_range(0.0..1.0); // Replace with actual precip-weighted sampling
    config.start_time + precip_weight * (config.end_time - config.start_time)
}

/// Precipitable water-weighted vertical sampling
pub fn parcel_release_height(config: &Config) -> f64 {
    let mut rng = rand::rng();
    let pw_weight = rng.random_range(0.0..1.0); // Replace with actual PW-weighted sampling
    config.constants.pressure_surface
        - pw_weight * (config.constants.pressure_surface - config.constants.pressure_top)
}

fn generate_random_parcels(config: &ReleaseConfig, count: usize) -> Vec<Parcel> {
    use rand::Rng;
    let mut rng = rand::rng();
    let mut parcels = Vec::new();

    for id in 0..count {
        let lon = rng.random_range(config.lon_min..=config.lon_max);
        let lat = rng.random_range(config.lat_min..=config.lat_max);
        let pressure = rng.random_range(config.pressure_min..=config.pressure_max);
        let time = rng.random_range(config.time_start..=config.time_end);

        parcels.push(Parcel::new(id as u32, lon, lat, pressure, time));
    }

    parcels
}

/// Generate single point release
fn generate_point_parcels(config: &ReleaseConfig) -> Vec<Parcel> {
    let lon = (config.lon_min + config.lon_max) / 2.0;
    let lat = (config.lat_min + config.lat_max) / 2.0;
    let pressure = (config.pressure_min + config.pressure_max) / 2.0;
    let time = config.time_start;

    vec![Parcel::new(0, lon, lat, pressure, time)]
}

/// Generate parcels along a line (e.g., flight track)
fn generate_line_parcels(config: &ReleaseConfig, waypoints: &[(f64, f64, f64)]) -> Vec<Parcel> {
    let mut parcels = Vec::new();
    let mut id = 0;

    for &(lon, lat, pressure) in waypoints {
        parcels.push(Parcel::new(id, lon, lat, pressure, config.time_start));
        id += 1;
    }

    parcels
}

/// Generate parcels uniformly over an area
fn generate_area_parcels(config: &ReleaseConfig) -> Vec<Parcel> {
    // For now, use grid pattern as area source
    generate_grid_parcels(config)
}
