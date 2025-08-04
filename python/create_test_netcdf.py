#!/usr/bin/env python3
"""
Create a test NetCDF file with WRF-like structure for testing the Rust NetCDF reader.
"""

import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
import os

def create_test_netcdf_file(filename="test_data/wrfout_d01_2023-07-31_12:00:00.nc"):
    """Create a test NetCDF file with WRF-like variables and dimensions."""
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(filename) if os.path.dirname(filename) else ".", exist_ok=True)
    
    # Dimensions
    nt = 3   # time steps
    nz = 5   # vertical levels (bottom_top)
    ny = 12  # south-north
    nx = 15  # west-east
    
    # Create NetCDF file
    with nc.Dataset(filename, 'w', format='NETCDF4') as ncfile:
        
        # Create dimensions
        ncfile.createDimension('Time', nt)
        ncfile.createDimension('bottom_top', nz)
        ncfile.createDimension('south_north', ny)
        ncfile.createDimension('west_east', nx)
        ncfile.createDimension('DateStrLen', 19)  # Create dimension for string length first
        
        # Create coordinate variables
        times = ncfile.createVariable('Times', 'S1', ('Time', 'DateStrLen'))
        times.description = 'times'
        
        # Create latitude and longitude
        xlat = ncfile.createVariable('XLAT', 'f4', ('south_north', 'west_east'))
        xlat.units = 'degree_north'
        xlat.description = 'Latitude'
        
        xlong = ncfile.createVariable('XLONG', 'f4', ('south_north', 'west_east'))
        xlong.units = 'degree_east'
        xlong.description = 'Longitude'
        
        # Create meteorological variables (4D: Time, bottom_top, south_north, west_east)
        u_var = ncfile.createVariable('U', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        u_var.units = 'm/s'
        u_var.description = 'U-component of wind'
        
        v_var = ncfile.createVariable('V', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        v_var.units = 'm/s'
        v_var.description = 'V-component of wind'
        
        w_var = ncfile.createVariable('W', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        w_var.units = 'm/s'
        w_var.description = 'W-component of wind'
        
        t_var = ncfile.createVariable('T', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        t_var.units = 'K'
        t_var.description = 'Temperature'
        
        qvapor_var = ncfile.createVariable('QVAPOR', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        qvapor_var.units = 'kg/kg'
        qvapor_var.description = 'Water vapor mixing ratio'
        
        rain_var = ncfile.createVariable('RAIN', 'f4', ('Time', 'south_north', 'west_east'), fill_value=-9999.0)
        rain_var.units = 'mm'
        rain_var.description = 'Accumulated rainfall'
        
        # Create pressure variable
        p_var = ncfile.createVariable('P', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        p_var.units = 'Pa'
        p_var.description = 'Pressure'
        
        # Fill with realistic test data
        
        # Create coordinate grids
        lons = np.linspace(-120.0, -110.0, nx)
        lats = np.linspace(30.0, 40.0, ny)
        lon_grid, lat_grid = np.meshgrid(lons, lats)
        
        xlat[:] = lat_grid
        xlong[:] = lon_grid
        
        # Create time strings
        base_time = datetime(2023, 7, 31, 12, 0, 0)
        time_strings = []
        for i in range(nt):
            time_str = (base_time + timedelta(hours=i)).strftime('%Y-%m-%d_%H:%M:%S')
            time_strings.append(time_str)
        
        # Fill time variable (convert strings to character arrays)
        for i, time_str in enumerate(time_strings):
            for j, char in enumerate(time_str):
                times[i, j] = char.encode('utf-8')
        
        # Create realistic meteorological data
        for t in range(nt):
            for k in range(nz):
                # Create some spatially and temporally varying patterns
                
                # Pressure decreases with height
                pressure_level = 100000.0 - k * 15000.0  # Pa
                p_var[t, k, :, :] = pressure_level * np.ones((ny, nx))
                
                # Temperature decreases with height and has spatial variation
                temp_base = 290.0 - k * 6.5  # K (lapse rate ~6.5 K/km)
                temp_variation = 5.0 * np.sin(2 * np.pi * lon_grid / 10.0) * np.cos(2 * np.pi * lat_grid / 10.0)
                t_var[t, k, :, :] = temp_base + temp_variation + np.random.normal(0, 0.5, (ny, nx))
                
                # U-wind component (westerly flow with variations)
                u_base = 10.0 + k * 2.0  # Stronger winds aloft
                u_variation = 3.0 * np.cos(2 * np.pi * lat_grid / 15.0)
                u_var[t, k, :, :] = u_base + u_variation + np.random.normal(0, 1.0, (ny, nx))
                
                # V-wind component (meridional flow)
                v_base = 2.0
                v_variation = 2.0 * np.sin(2 * np.pi * lon_grid / 20.0)
                v_var[t, k, :, :] = v_base + v_variation + np.random.normal(0, 0.8, (ny, nx))
                
                # W-wind component (small vertical motions)
                w_var[t, k, :, :] = 0.1 * np.sin(2 * np.pi * lon_grid / 8.0) * np.sin(2 * np.pi * lat_grid / 8.0) + np.random.normal(0, 0.05, (ny, nx))
                
                # Water vapor mixing ratio (decreases with height)
                qvapor_base = 0.015 * np.exp(-k * 0.5)  # kg/kg
                qvapor_var[t, k, :, :] = qvapor_base * np.ones((ny, nx)) + np.random.normal(0, 0.001, (ny, nx))
            
            # Rainfall (2D field)
            rain_pattern = 2.0 * np.exp(-((lon_grid + 115.0)**2 + (lat_grid - 35.0)**2) / 10.0)
            rain_var[t, :, :] = rain_pattern + np.random.exponential(0.5, (ny, nx))
        
        # Add global attributes
        ncfile.title = "Test WRF-like NetCDF file for Rust reader testing"
        ncfile.institution = "Test Institution"
        ncfile.source = "create_test_netcdf.py"
        ncfile.history = f"Created on {datetime.now().isoformat()}"
        ncfile.Conventions = "CF-1.6"
        
        print(f"Created test NetCDF file: {filename}")
        print(f"Dimensions: Time={nt}, bottom_top={nz}, south_north={ny}, west_east={nx}")
        print(f"Variables: U, V, W, T, QVAPOR, RAIN, P, XLAT, XLONG")

def create_smaller_test_file(filename="test_data/small_wrfout_d01_2023-07-31_12:00:00.nc"):
    """Create a smaller test file for boundary trimming tests."""
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(filename) if os.path.dirname(filename) else ".", exist_ok=True)
    
    # Small dimensions for testing boundary trimming
    nt = 2   # time steps
    nz = 3   # vertical levels
    ny = 6   # south-north (will be trimmed to 4)
    nx = 8   # west-east (will be trimmed to 6)
    
    with nc.Dataset(filename, 'w', format='NETCDF4') as ncfile:
        
        # Create dimensions
        ncfile.createDimension('Time', nt)
        ncfile.createDimension('bottom_top', nz)
        ncfile.createDimension('south_north', ny)
        ncfile.createDimension('west_east', nx)
        
        # Create a simple U variable for testing
        u_var = ncfile.createVariable('U', 'f4', ('Time', 'bottom_top', 'south_north', 'west_east'), fill_value=-9999.0)
        u_var.units = 'm/s'
        u_var.description = 'U-component of wind'
        
        # Fill with simple test pattern
        for t in range(nt):
            for k in range(nz):
                u_var[t, k, :, :] = np.arange(ny * nx).reshape(ny, nx) + t * 100 + k * 10
        
        print(f"Created small test NetCDF file: {filename}")
        print(f"Dimensions: Time={nt}, bottom_top={nz}, south_north={ny}, west_east={nx}")

if __name__ == "__main__":
    create_test_netcdf_file()
    create_smaller_test_file()
    print("\nTest NetCDF files created successfully!")
    print("You can now run the Rust tests with actual NetCDF data.")
