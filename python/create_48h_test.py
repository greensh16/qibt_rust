#!/usr/bin/env python3
"""
Create a test NetCDF file with 48 hours of 10-minute data for testing SingleFileDataLoader
"""

import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta

def create_48h_test_netcdf(filename='test_data/single_file_test.nc'):
    """Create a NetCDF file with 48 hours of 10-minute interval data"""
    
    # Create NetCDF file
    with nc.Dataset(filename, 'w', format='NETCDF4') as ds:
        # Create dimensions
        time_dim = ds.createDimension('Time', None)  # unlimited - WRF style
        level_dim = ds.createDimension('bottom_top', 10)  # WRF style
        lat_dim = ds.createDimension('south_north', 20)   # WRF style
        lon_dim = ds.createDimension('west_east', 25)     # WRF style
        
        # Create coordinate variables
        times = ds.createVariable('XTIME', 'f8', ('Time',))
        levels = ds.createVariable('bottom_top', 'i4', ('bottom_top',))
        lats = ds.createVariable('XLAT', 'f4', ('south_north', 'west_east'))
        lons = ds.createVariable('XLONG', 'f4', ('south_north', 'west_east'))
        
        # Create meteorological variables (WRF style dimensions: Time, bottom_top, south_north, west_east)
        temperature = ds.createVariable('T', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        u_wind = ds.createVariable('U', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        v_wind = ds.createVariable('V', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        w_wind = ds.createVariable('W', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        qvapor = ds.createVariable('QVAPOR', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        
        # Generate coordinate data
        lat_1d = np.linspace(30.0, 40.0, 20)  # Realistic latitude range
        lon_1d = np.linspace(-120.0, -110.0, 25)  # Realistic longitude range
        lat_2d, lon_2d = np.meshgrid(lat_1d, lon_1d, indexing='ij')
        level_data = np.arange(1, 11)
        
        # Set coordinate data
        lats[:] = lat_2d
        lons[:] = lon_2d
        levels[:] = level_data
        
        # Generate 48 hours of 10-minute data (288 timesteps)
        start_time = datetime(2025, 8, 1, 0, 0, 0)
        end_time = start_time + timedelta(hours=48)
        
        time_list = []
        current_time = start_time
        
        print(f"Generating time series from {start_time} to {end_time}...")
        
        while current_time < end_time:
            # XTIME in WRF is typically minutes since simulation start
            minutes_since_start = (current_time - start_time).total_seconds() / 60.0
            time_list.append(minutes_since_start)
            current_time += timedelta(minutes=10)
        
        # Set time data
        times[:] = time_list
        num_times = len(time_list)
        
        print(f"Generated {num_times} timesteps...")
        print("Creating synthetic meteorological data...")
        
        # Generate synthetic meteorological data with realistic values and variation
        np.random.seed(42)  # For reproducible results
        
        # Temperature: varies between 280-310K with diurnal cycle
        temp_base = np.random.uniform(280, 310, size=(num_times, 10, 20, 25))
        for t in range(num_times):
            # Add diurnal cycle (24-hour period)
            diurnal_variation = 10 * np.sin(2 * np.pi * t / 144)  # 144 timesteps = 24 hours
            temp_base[t] += diurnal_variation
            # Add altitude variation (temperature decreases with height)
            for k in range(10):
                temp_base[t, k] -= k * 5  # 5K decrease per level
        
        temperature[:] = temp_base
        
        # Wind components: realistic wind speeds with some structure
        u_base = np.random.uniform(-15, 15, size=(num_times, 10, 20, 25))
        v_base = np.random.uniform(-15, 15, size=(num_times, 10, 20, 25))
        w_base = np.random.uniform(-2, 2, size=(num_times, 10, 20, 25))  # Smaller vertical velocities
        
        # Add some coherent wind patterns
        for t in range(num_times):
            for k in range(10):
                # Add some spatial correlation
                u_base[t, k] += 5 * np.sin(lat_2d * np.pi / 180) * np.cos(lon_2d * np.pi / 180)
                v_base[t, k] += 3 * np.cos(lat_2d * np.pi / 180) * np.sin(lon_2d * np.pi / 180)
        
        u_wind[:] = u_base
        v_wind[:] = v_base
        w_wind[:] = w_base
        
        # Water vapor mixing ratio: realistic values 0-20 g/kg
        qvapor_base = np.random.uniform(0.001, 0.020, size=(num_times, 10, 20, 25))
        # Decrease with height
        for k in range(10):
            qvapor_base[:, k] *= np.exp(-k * 0.3)
        qvapor[:] = qvapor_base
        
        # Set attributes
        ds.title = '48-hour test NetCDF file for SingleFileDataLoader'
        ds.institution = 'Test Data Generator'
        ds.source = 'Synthetic meteorological data'
        ds.history = f'Created on {datetime.now().isoformat()}'
        ds.Conventions = 'CF-1.6'
        
        # Variable attributes
        times.long_name = 'minutes since simulation start'
        times.units = 'minutes'
        
        temperature.long_name = 'perturbation potential temperature (theta-t0)'
        temperature.units = 'K'
        temperature.description = 'POTENTIAL TEMPERATURE'
        
        u_wind.long_name = 'x-wind component'
        u_wind.units = 'm s-1'
        u_wind.description = 'U'
        
        v_wind.long_name = 'y-wind component'
        v_wind.units = 'm s-1'
        v_wind.description = 'V'
        
        w_wind.long_name = 'z-wind component'
        w_wind.units = 'm s-1'
        w_wind.description = 'W'
        
        qvapor.long_name = 'Water vapor mixing ratio'
        qvapor.units = 'kg kg-1'
        qvapor.description = 'QVAPOR'
        
        lats.long_name = 'LATITUDE, SOUTH IS NEGATIVE'
        lats.units = 'degree_north'
        lons.long_name = 'LONGITUDE, WEST IS NEGATIVE'
        lons.units = 'degree_east'
        
    print(f"Created test NetCDF file: {filename}")
    print(f"Dimensions: Time={num_times}, bottom_top=10, south_north=20, west_east=25")
    print(f"Time range: {start_time} to {end_time}")
    print(f"Time step: 10 minutes")
    print(f"Total duration: 48 hours")
    print(f"File size: {filename}")

if __name__ == '__main__':
    # Ensure test_data directory exists
    import os
    os.makedirs('test_data', exist_ok=True)
    
    create_48h_test_netcdf()
