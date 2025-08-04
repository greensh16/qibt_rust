pub fn julian(year: i64, month: i64, day: i64) -> i64 {
    day - 32075
        + 1461 * (year + 4800 + (month - 14) / 12) / 4
        + 367 * (month - 2 - (month - 14) / 12 * 12) / 12
        - 3 * ((year + 4900 + (month - 14) / 12) / 100) / 4
}

pub fn gregorian(jd: i64) -> (i64, i64, i64) {
    let l = jd + 68569;
    let n = 4 * l / 146097;
    let l = l - (146097 * n + 3) / 4;
    let i = 4000 * (l + 1) / 1461001;
    let l = l - 1461 * i / 4 + 31;
    let j = 80 * l / 2447;
    let k = l - 2447 * j / 80;
    let l = j / 11;
    let j = j + 2 - 12 * l;
    let i = 100 * (n - 49) + i + l;

    (i, j, k)
}

pub fn simlength(
    start_day: i64,
    start_month: i64,
    start_year: i64,
    end_day: i64,
    end_month: i64,
    end_year: i64,
) -> i64 {
    let start_jd = julian(start_year, start_month, start_day);
    let end_jd = julian(end_year, end_month, end_day);
    end_jd - start_jd
}

pub fn day_month_year(
    sim_day: i64,
    start_day: i64,
    start_month: i64,
    start_year: i64,
) -> (i64, i64, i64) {
    if sim_day == 1 {
        (start_year, start_month, start_day)
    } else {
        let jday = julian(start_year, start_month, start_day) + (sim_day - 1);
        gregorian(jday)
    }
}

use std::fmt;

/// Gregorian date representation
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GregorianDate {
    pub year: i32,
    pub month: u8,
    pub day: u8,
    pub hour: u8,
    pub minute: u8,
    pub second: u8,
}

impl fmt::Display for GregorianDate {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:04}-{:02}-{:02} {:02}:{:02}:{:02}",
            self.year, self.month, self.day, self.hour, self.minute, self.second
        )
    }
}

/// Convert Julian day number to Gregorian date
pub fn julian_to_gregorian(julian_day: f64) -> GregorianDate {
    // Split into integer day and fractional part
    let jd_int = julian_day.floor() as i64;
    let jd_frac = julian_day - jd_int as f64;

    // Convert to Gregorian calendar using the algorithm from
    // "Astronomical Algorithms" by Jean Meeus
    let a = jd_int + 32044;
    let b = (4 * a + 3) / 146097;
    let c = a - (146097 * b) / 4;
    let d = (4 * c + 3) / 1461;
    let e = c - (1461 * d) / 4;
    let m = (5 * e + 2) / 153;

    let day = e - (153 * m + 2) / 5 + 1;
    let month = m + 3 - 12 * (m / 10);
    let year = 100 * b + d - 4800 + m / 10;

    // Convert fractional day to hours, minutes, seconds
    let total_seconds = (jd_frac * 86400.0).round() as u32;
    let hour = (total_seconds / 3600) as u8;
    let minute = ((total_seconds % 3600) / 60) as u8;
    let second = (total_seconds % 60) as u8;

    GregorianDate {
        year: year as i32,
        month: month as u8,
        day: day as u8,
        hour,
        minute,
        second,
    }
}

/// Convert Gregorian date to Julian day number
pub fn gregorian_to_julian(date: &GregorianDate) -> f64 {
    let y = if date.month <= 2 {
        date.year - 1
    } else {
        date.year
    };

    let m = if date.month <= 2 {
        date.month as i32 + 12
    } else {
        date.month as i32
    };

    let a = y / 100;
    let b = 2 - a + a / 4; // Gregorian calendar correction

    let jd_int = (365.25 * (y + 4716) as f64).floor() as i64
        + (30.6001 * (m + 1) as f64).floor() as i64
        + date.day as i64
        + b as i64
        - 1524;

    // Add fractional day from time
    let day_fraction =
        (date.hour as f64 * 3600.0 + date.minute as f64 * 60.0 + date.second as f64) / 86400.0;

    jd_int as f64 + day_fraction
}

/// Calculate the number of days in a given month
pub fn days_in_month(year: i32, month: u8) -> u8 {
    match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 => {
            if is_leap_year(year) {
                29
            } else {
                28
            }
        }
        _ => panic!("Invalid month: {}", month),
    }
}

/// Check if a year is a leap year
pub fn is_leap_year(year: i32) -> bool {
    (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)
}

/// Add time interval to Julian day
pub fn add_time_to_julian(julian_day: f64, hours: f64) -> f64 {
    julian_day + hours / 24.0
}

/// Calculate time difference in hours between two Julian days
pub fn julian_time_diff_hours(jd1: f64, jd2: f64) -> f64 {
    (jd2 - jd1) * 24.0
}

/// Parse time string in various formats
pub fn parse_time_string(time_str: &str) -> Result<GregorianDate, String> {
    // Try different formats
    if let Ok(date) = parse_iso_format(time_str) {
        return Ok(date);
    }

    if let Ok(date) = parse_simple_format(time_str) {
        return Ok(date);
    }

    Err(format!("Could not parse time string: {}", time_str))
}

/// Parse ISO 8601 format: YYYY-MM-DDTHH:MM:SS or YYYY-MM-DD HH:MM:SS
fn parse_iso_format(time_str: &str) -> Result<GregorianDate, String> {
    let parts: Vec<&str> = time_str.split(&['T', ' '][..]).collect();
    if parts.len() != 2 {
        return Err("Invalid ISO format".to_string());
    }

    let date_parts: Vec<&str> = parts[0].split('-').collect();
    let time_parts: Vec<&str> = parts[1].split(':').collect();

    if date_parts.len() != 3 || time_parts.len() != 3 {
        return Err("Invalid date/time components".to_string());
    }

    Ok(GregorianDate {
        year: date_parts[0].parse().map_err(|_| "Invalid year")?,
        month: date_parts[1].parse().map_err(|_| "Invalid month")?,
        day: date_parts[2].parse().map_err(|_| "Invalid day")?,
        hour: time_parts[0].parse().map_err(|_| "Invalid hour")?,
        minute: time_parts[1].parse().map_err(|_| "Invalid minute")?,
        second: time_parts[2].parse().map_err(|_| "Invalid second")?,
    })
}

/// Parse simple format: YYYYMMDDHH
fn parse_simple_format(time_str: &str) -> Result<GregorianDate, String> {
    if time_str.len() != 10 {
        return Err("Simple format must be YYYYMMDDHH".to_string());
    }

    Ok(GregorianDate {
        year: time_str[0..4].parse().map_err(|_| "Invalid year")?,
        month: time_str[4..6].parse().map_err(|_| "Invalid month")?,
        day: time_str[6..8].parse().map_err(|_| "Invalid day")?,
        hour: time_str[8..10].parse().map_err(|_| "Invalid hour")?,
        minute: 0,
        second: 0,
    })
}

/// Calculate simulation end time given start time and duration
pub fn calculate_end_time(start_julian: f64, duration_hours: f64) -> f64 {
    start_julian + duration_hours / 24.0
}

/// Generate time sequence for trajectory integration
pub fn generate_time_sequence(
    start_julian: f64,
    end_julian: f64,
    time_step_hours: f64,
) -> Vec<f64> {
    let mut times = Vec::new();
    let time_step_days = time_step_hours / 24.0;

    let mut current_time = start_julian;

    // For backward trajectories, we integrate backward in time
    if end_julian < start_julian {
        while current_time >= end_julian {
            times.push(current_time);
            current_time -= time_step_days;
        }
    } else {
        // Forward trajectories
        while current_time <= end_julian {
            times.push(current_time);
            current_time += time_step_days;
        }
    }

    times
}

/// Format Julian day as string
pub fn format_julian_day(julian_day: f64) -> String {
    let date = julian_to_gregorian(julian_day);
    format!("{}", date)
}

/// Get current Julian day (system time)
pub fn current_julian_day() -> f64 {
    use std::time::{SystemTime, UNIX_EPOCH};

    let unix_time = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards")
        .as_secs() as f64;

    // Unix epoch (1970-01-01 00:00:00) corresponds to Julian day 2440587.5
    2440587.5 + unix_time / 86400.0
}

/// Time zone conversion (placeholder - assumes UTC)
pub fn convert_timezone(julian_day: f64, _timezone_offset_hours: f64) -> f64 {
    // TODO: Implement proper timezone conversion
    // For now, just return the input (assumes all times are UTC)
    julian_day
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_julian_gregorian_round_trip() {
        for year in 1900..2100 {
            for month in 1..=12 {
                for day in 1..=days_in_month(year, month) {
                    let jd = julian(year as i64, month as i64, day as i64);
                    let (y, m, d) = gregorian(jd);
                    assert_eq!((year as i64, month as i64, day as i64), (y, m, d));
                }
            }
        }
    }

    #[test]
    fn test_julian_and_gregorian_specific() {
        let date = GregorianDate {
            year: 2000,
            month: 1,
            day: 1,
            hour: 0,
            minute: 0,
            second: 0,
        };

        let jd = gregorian_to_julian(&date);
        let back_to_gregorian = julian_to_gregorian(jd);

        assert_eq!(date, back_to_gregorian);
    }

    #[test]
    fn test_simlength() {
        let start_day = 1;
        let start_month = 1;
        let start_year = 2023;
        let end_day = 31;
        let end_month = 12;
        let end_year = 2023;

        assert_eq!(
            simlength(
                start_day,
                start_month,
                start_year,
                end_day,
                end_month,
                end_year
            ),
            364
        );
    }

    #[test]
    fn test_day_month_year() {
        let start_day = 1;
        let start_month = 1;
        let start_year = 2023;

        let sim_day = 365;
        let (year, month, day) = day_month_year(sim_day, start_day, start_month, start_year);

        assert_eq!((year, month, day), (2023, 12, 31));
    }

    #[test]
    fn test_leap_year() {
        assert!(is_leap_year(2000));
        assert!(!is_leap_year(1900));
        assert!(is_leap_year(2004));
        assert!(!is_leap_year(2001));
    }
}
