#!/usr/bin/env python3
"""
Create test Zarr and NetCDF files with identical data for parity testing.
"""

import numpy as np
import zarr
import netCDF4 as nc
import os
import shutil
from pathlib import Path

def create_test_data():
    """Create consistent test data for both formats."""
    # Dimensions
    nt, nz, ny, nx = 3, 5, 10, 12
    
    # Coordinate data
    times = np.arange(nt, dtype=np.float32)
    levels = np.array([1000.0, 850.0, 700.0, 500.0, 300.0], dtype=np.float32)
    
    # Create 2D coordinate grids
    lats = np.linspace(30.0, 40.0, ny, dtype=np.float32)
    lons = np.linspace(-120.0, -110.0, nx, dtype=np.float32)
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    
    # Create 4D meteorological variables
    shape_4d = (nt, nz, ny, nx)
    total_elements = nt * nz * ny * nx
    
    # U-component wind - linear pattern for easy verification
    u_data = np.arange(total_elements, dtype=np.float32).reshape(shape_4d) * 0.1 + 10.0
    
    # V-component wind
    v_data = np.arange(total_elements, dtype=np.float32).reshape(shape_4d) * 0.05 + 5.0
    
    # Temperature
    temp_data = np.arange(total_elements, dtype=np.float32).reshape(shape_4d) * 0.01 + 280.0
    
    return {
        'dimensions': {'time': nt, 'level': nz, 'lat': ny, 'lon': nx},
        'coordinates': {
            'time': times,
            'level': levels,
            'lat': lat_grid,
            'lon': lon_grid
        },
        'variables': {
            'u': u_data,
            'v': v_data,
            'temp': temp_data
        }
    }

def create_zarr_store(zarr_path, data):
    """Create a Zarr store with test data."""
    if os.path.exists(zarr_path):
        shutil.rmtree(zarr_path)
    
    # Create root group
    root = zarr.open_group(zarr_path, mode='w')
    
    # Add global attributes
    root.attrs['title'] = 'Test Zarr store for parity testing'
    root.attrs['institution'] = 'Test Institution'
    root.attrs['source'] = 'create_test_zarr.py'
    root.attrs['Conventions'] = 'CF-1.6'
    
    # Create coordinate variables
    time_var = root.create_dataset('time', data=data['coordinates']['time'],
                                   chunks=True, compressor=zarr.Blosc())
    time_var.attrs['units'] = 'hours since 2023-01-01'
    time_var.attrs['long_name'] = 'Time'
    
    level_var = root.create_dataset('level', data=data['coordinates']['level'],
                                    chunks=True, compressor=zarr.Blosc())
    level_var.attrs['units'] = 'hPa'
    level_var.attrs['long_name'] = 'Pressure level'
    
    lat_var = root.create_dataset('lat', data=data['coordinates']['lat'],
                                  chunks=True, compressor=zarr.Blosc())
    lat_var.attrs['units'] = 'degree_north'
    lat_var.attrs['long_name'] = 'Latitude'
    
    lon_var = root.create_dataset('lon', data=data['coordinates']['lon'],
                                  chunks=True, compressor=zarr.Blosc())
    lon_var.attrs['units'] = 'degree_east'
    lon_var.attrs['long_name'] = 'Longitude'
    
    # Create meteorological variables
    u_var = root.create_dataset('u', data=data['variables']['u'],
                                chunks=(1, 2, 5, 6), compressor=zarr.Blosc())
    u_var.attrs['units'] = 'm/s'
    u_var.attrs['long_name'] = 'U-component of wind'
    
    v_var = root.create_dataset('v', data=data['variables']['v'],
                                chunks=(1, 2, 5, 6), compressor=zarr.Blosc())
    v_var.attrs['units'] = 'm/s'
    v_var.attrs['long_name'] = 'V-component of wind'
    
    temp_var = root.create_dataset('temp', data=data['variables']['temp'],
                                   chunks=(1, 2, 5, 6), compressor=zarr.Blosc())
    temp_var.attrs['units'] = 'K'
    temp_var.attrs['long_name'] = 'Temperature'
    
    print(f"Created Zarr store: {zarr_path}")
    print(f"  Dimensions: {data['dimensions']}")
    print(f"  Variables: {list(data['variables'].keys())}")

def create_netcdf_file(netcdf_path, data):
    """Create a NetCDF file with test data."""
    if os.path.exists(netcdf_path):
        os.remove(netcdf_path)
    
    with nc.Dataset(netcdf_path, 'w', format='NETCDF4') as ncfile:
        # Create dimensions
        ncfile.createDimension('time', data['dimensions']['time'])
        ncfile.createDimension('level', data['dimensions']['level'])
        ncfile.createDimension('lat', data['dimensions']['lat'])
        ncfile.createDimension('lon', data['dimensions']['lon'])
        
        # Add global attributes
        ncfile.title = 'Test NetCDF file for parity testing'
        ncfile.institution = 'Test Institution'
        ncfile.source = 'create_test_zarr.py'
        ncfile.Conventions = 'CF-1.6'
        
        # Create coordinate variables
        time_var = ncfile.createVariable('time', 'f4', ('time',))
        time_var.units = 'hours since 2023-01-01'
        time_var.long_name = 'Time'
        time_var[:] = data['coordinates']['time']
        
        level_var = ncfile.createVariable('level', 'f4', ('level',))
        level_var.units = 'hPa'
        level_var.long_name = 'Pressure level'
        level_var[:] = data['coordinates']['level']
        
        # Create coordinate grids (NetCDF style with WRF names)
        xlat_var = ncfile.createVariable('XLAT', 'f4', ('lat', 'lon'))
        xlat_var.units = 'degree_north'
        xlat_var.long_name = 'Latitude'
        xlat_var[:] = data['coordinates']['lat']
        
        xlong_var = ncfile.createVariable('XLONG', 'f4', ('lat', 'lon'))
        xlong_var.units = 'degree_east'
        xlong_var.long_name = 'Longitude'
        xlong_var[:] = data['coordinates']['lon']
        
        # Create meteorological variables (using WRF-style names)
        u_var = ncfile.createVariable('U', 'f4', ('time', 'level', 'lat', 'lon'),
                                      fill_value=-9999.0)
        u_var.units = 'm/s'
        u_var.long_name = 'U-component of wind'
        u_var[:] = data['variables']['u']
        
        v_var = ncfile.createVariable('V', 'f4', ('time', 'level', 'lat', 'lon'),
                                      fill_value=-9999.0)
        v_var.units = 'm/s'
        v_var.long_name = 'V-component of wind'
        v_var[:] = data['variables']['v']
        
        t_var = ncfile.createVariable('T', 'f4', ('time', 'level', 'lat', 'lon'),
                                      fill_value=-9999.0)
        t_var.units = 'K'
        t_var.long_name = 'Temperature'
        t_var[:] = data['variables']['temp']
    
    print(f"Created NetCDF file: {netcdf_path}")
    print(f"  Dimensions: {data['dimensions']}")
    print(f"  Variables: U, V, T, XLAT, XLONG")

def verify_parity(zarr_path, netcdf_path):
    """Verify that Zarr and NetCDF files contain identical data."""
    print("\nVerifying data parity...")
    
    # Open Zarr store
    zarr_root = zarr.open(zarr_path, mode='r')
    
    # Open NetCDF file
    with nc.Dataset(netcdf_path, 'r') as ncfile:
        # Compare U variable (zarr 'u' vs netcdf 'U')
        zarr_u = zarr_root['u'][:]
        netcdf_u = ncfile.variables['U'][:]
        
        if np.allclose(zarr_u, netcdf_u):
            print("✓ U/u variable data matches")
        else:
            print("✗ U/u variable data differs")
            print(f"  Max difference: {np.max(np.abs(zarr_u - netcdf_u))}")
        
        # Compare coordinate data
        zarr_lat = zarr_root['lat'][:]
        netcdf_lat = ncfile.variables['XLAT'][:]
        
        if np.allclose(zarr_lat, netcdf_lat):
            print("✓ Latitude data matches")
        else:
            print("✗ Latitude data differs")
        
        # Check dimensions
        zarr_shape = zarr_u.shape
        netcdf_shape = netcdf_u.shape
        
        if zarr_shape == netcdf_shape:
            print(f"✓ Array shapes match: {zarr_shape}")
        else:
            print(f"✗ Array shapes differ: Zarr {zarr_shape} vs NetCDF {netcdf_shape}")

def main():
    """Main function to create test data files."""
    # Create output directories
    test_data_dir = Path("test_data")
    test_data_dir.mkdir(exist_ok=True)
    
    zarr_path = test_data_dir / "parity_test.zarr"
    netcdf_path = test_data_dir / "parity_test.nc"
    
    # Generate consistent test data
    print("Generating test data...")
    data = create_test_data()
    
    # Create both formats
    create_zarr_store(str(zarr_path), data)
    create_netcdf_file(str(netcdf_path), data)
    
    # Verify parity
    verify_parity(str(zarr_path), str(netcdf_path))
    
    print(f"\nTest files created successfully!")
    print(f"Zarr store: {zarr_path}")
    print(f"NetCDF file: {netcdf_path}")
    print("\nThese files contain identical data and can be used for parity testing.")

if __name__ == "__main__":
    main()
