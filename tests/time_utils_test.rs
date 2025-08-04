use qibt_rust::time_utils::{
    day_month_year, days_in_month, gregorian, is_leap_year, julian, simlength,
};

#[test]
fn test_julian_gregorian_round_trip() {
    // Test a subset of dates to avoid long test times
    for year in [1900, 1950, 2000, 2023, 2024] {
        for month in 1..=12 {
            for day in 1..=days_in_month(year, month) {
                let jd = julian(year as i64, month as i64, day as i64);
                let (y, m, d) = gregorian(jd);
                assert_eq!(
                    (year as i64, month as i64, day as i64),
                    (y, m, d),
                    "Failed for {}-{}-{}",
                    year,
                    month,
                    day
                );
            }
        }
    }
}

#[test]
fn test_specific_julian_dates() {
    // Test Y2K date
    let jd_y2k = julian(2000, 1, 1);
    let (y, m, d) = gregorian(jd_y2k);
    assert_eq!((y, m, d), (2000, 1, 1));

    // Test epoch date (1970-01-01)
    let jd_epoch = julian(1970, 1, 1);
    let (y, m, d) = gregorian(jd_epoch);
    assert_eq!((y, m, d), (1970, 1, 1));

    // Test leap year date
    let jd_leap = julian(2024, 2, 29);
    let (y, m, d) = gregorian(jd_leap);
    assert_eq!((y, m, d), (2024, 2, 29));
}

#[test]
fn test_simlength() {
    // Test one year simulation
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

    // Test leap year
    let leap_year_days = simlength(1, 1, 2024, 31, 12, 2024);
    assert_eq!(leap_year_days, 365); // 2024 is a leap year
}

#[test]
fn test_day_month_year() {
    let start_day = 1;
    let start_month = 1;
    let start_year = 2023;

    // First day should return start date
    let (year, month, day) = day_month_year(1, start_day, start_month, start_year);
    assert_eq!((year, month, day), (2023, 1, 1));

    // Last day of year
    let (year, month, day) = day_month_year(365, start_day, start_month, start_year);
    assert_eq!((year, month, day), (2023, 12, 31));

    // Test day 32 (February 1st)
    let (year, month, day) = day_month_year(32, start_day, start_month, start_year);
    assert_eq!((year, month, day), (2023, 2, 1));
}

#[test]
fn test_leap_year() {
    assert!(is_leap_year(2000)); // Divisible by 400
    assert!(!is_leap_year(1900)); // Divisible by 100, not by 400
    assert!(is_leap_year(2004)); // Divisible by 4, not by 100
    assert!(!is_leap_year(2001)); // Not divisible by 4
    assert!(is_leap_year(2024)); // Current leap year
}

#[test]
fn test_days_in_month() {
    // Test regular year
    assert_eq!(days_in_month(2023, 1), 31); // January
    assert_eq!(days_in_month(2023, 2), 28); // February (non-leap)
    assert_eq!(days_in_month(2023, 4), 30); // April
    assert_eq!(days_in_month(2023, 12), 31); // December

    // Test leap year
    assert_eq!(days_in_month(2024, 2), 29); // February (leap year)
}

#[test]
fn test_month_boundaries() {
    // Test crossing month boundaries
    let jd_jan31 = julian(2023, 1, 31);
    let jd_feb1 = julian(2023, 2, 1);
    assert_eq!(jd_feb1 - jd_jan31, 1);

    // Test crossing year boundaries
    let jd_dec31 = julian(2023, 12, 31);
    let jd_jan1 = julian(2024, 1, 1);
    assert_eq!(jd_jan1 - jd_dec31, 1);
}
