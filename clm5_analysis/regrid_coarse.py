# -*- coding: utf-8 -*-
"""
Desc: for revisions we want to analyse everything at 0.25 and/or 0.5 degrees -> upscale 1km simulations
-> use xesmf_env_new2 !!!
Created on 29.01.24 08:54
@author: malle
"""
import pandas as pd
from pathlib import Path
import xarray as xr
import glob
import xesmf as xe
import scipy.io
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import platform


if platform.system() == 'Linux':
    bf = Path('/home/lud11/malle/CLM5_CH')
else:
    bf = Path('L:\malle\CLM5_CH')


bf_jim = bf / 'FSM_new' / 'analysed_grid_python'
bf_jim_out = Path('/home/lud11/malle/CLM5_CH/FSM_new/analysed_grid_python_out')
# import surfdata & create target grid for xesmf:
surf_in = xr.open_dataset(bf / 'surfdata_1km_CH_v3_hist_16pfts_Irrig_CMIP6_NEW.nc')
lat_new_1km = surf_in.LATIXY.data[:, 0]
lon_new_1km = surf_in.LONGXY.data[0, :]
ds_target_1km = xr.Dataset({"lat": (["lat"], lat_new_1km), "lon": (["lon"], lon_new_1km)})

surf_in = xr.open_dataset(bf / 'surfdata_CH_025deg_v2_hist_16pfts_Irrig_CMIP6_simyr2000_c230113.nc')
lat_new = surf_in.LATIXY.data[:, 0]
lon_new = surf_in.LONGXY.data[0, :]
ds_target_025 = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

surf_in = xr.open_dataset(bf / 'surfdata_CH_05deg_v2_hist_16pfts_Irrig_CMIP6_simyr2000_c230113.nc')
lat_new = surf_in.LATIXY.data[:, 0]
lon_new = surf_in.LONGXY.data[0, :]
ds_target_05 = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

regridder_05 = xe.Regridder(ds_target_1km, ds_target_05, "conservative")
regridder_025 = xe.Regridder(ds_target_1km, ds_target_025, "conservative")
regridder_coarse_025 = xe.Regridder(ds_target_05, ds_target_025, "nearest_s2d")


# loop through 4 seasons of interest
for file_jim in list(bf_jim.iterdir()):
    mat = scipy.io.loadmat(file_jim, simplify_cells=True)
    mat_jim = mat["OSHD"]
    time_stamp = pd.to_datetime(mat_jim["time_stamps"])
    hs_jim = mat_jim['hs']

    da = xr.DataArray(
        data=np.flipud(hs_jim),
        coords=[("lat", lat_new_1km.astype("float32")),
                ("lon", lon_new_1km.astype("float32")),
                ("time", time_stamp)]
        )
    ds_out = xr.Dataset(data_vars={"DATA": da})

    ds_out.to_netcdf(f"{bf_jim_out}/analysis_{str(time_stamp[-1].year)}_1km.nc")

    snow_oshd_025 = regridder_025(ds_out)
    snow_oshd_025.to_netcdf(f"{bf_jim_out}/analysis_{str(time_stamp[-1].year)}_025deg.nc")



run_in_all = ['OSHD_FILES', 'OSHD_FILES_OLD', 'CRUJRA_FILES_noLapse', 'CRUJRA_FILES_noLapse_OLD', 'CRUJRA_FILES',
              'CRUJRA_FILES_OLD']

run_in_05 = ['OSHD_FILES_05_new', 'CRUJRA_FILES_05deg_cru_new', 'CRUJRA_FILES_05deg_cru_new_lapse']

# now loop through all 05 simulations... -> to 0.25degree
for run_in in run_in_05:
    print(run_in)
    snow_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/SNOW_DEPTH*')[0])
    evap_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/QFLX_EVAP_TOT*')[0]) * 3600 * 24  # mm/s to mm/day
    evap_in_sum_yr = evap_in.resample(time='1Y').sum()

    snow_out_025 = regridder_coarse_025(snow_in)
    file_snow_025 = bf / run_in / 'downscale_025_SNOW_DEPTH.nc'
    snow_out_025.to_netcdf(file_snow_025)

    # now evapo
    evap_out_05 = regridder_coarse_025(evap_in_sum_yr)
    file_05 = bf / run_in / 'downscale_05_QFLX_EVAP_TOT_YR_SUM.nc'
    evap_out_05.to_netcdf(file_05)


# now loop through all 1km simulations... -> to 0.25degree
for run_in in run_in_all:
    print(run_in)
    snow_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/SNOW_DEPTH*')[0])
    snow_in = snow_in.where(~np.isnan(da.isel(time=0)), np.nan)  # set to nan where OSHD is set to 0 as well!
    evap_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/QFLX_EVAP_TOT*')[0]) * 3600 * 24  # mm/s to mm/day
    evap_in_sum_yr = evap_in.resample(time='1Y').sum()

    snow_out_025 = regridder_025(snow_in)
    file_snow_025 = bf / run_in / 'upscale_025_SNOW_DEPTH.nc'
    snow_out_025.to_netcdf(file_snow_025)

    snow_out_05 = regridder_05(snow_in)
    file_snow_05 = bf / run_in / 'upscale_05_SNOW_DEPTH.nc'
    snow_out_05.to_netcdf(file_snow_05)

    # now 0.25degree ET
    evap_out_025 = regridder_025(evap_in_sum_yr)
    file_025 = bf / run_in / 'upscale_025_QFLX_EVAP_TOT_YR_SUM.nc'
    evap_out_025.to_netcdf(file_025)

    # now 0.5degree
    evap_out_05 = regridder_05(evap_in_sum_yr)
    file_05 = bf / run_in / 'upscale_05_QFLX_EVAP_TOT_YR_SUM.nc'
    evap_out_05.to_netcdf(file_05)
