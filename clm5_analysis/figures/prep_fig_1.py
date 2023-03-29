# -*- coding: utf-8 -*-
"""
Desc:
Created on 08.12.22 18:02
@author: malle
"""

import numpy as np
from pathlib import Path
import xarray as xr
import matplotlib
import platform
import os

# set font type for PDFs and EPSs
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# make switch to windows if working from home
if platform.system() == 'Windows':
    base_dir = Path('L:\malle\CLM5_CH')
else:
    base_dir = Path('/home/lud11/malle/CLM5_CH')
    import xesmf as xe

# load domain, surface data (for vegetation cover and regridding) and temperatures
domain = xr.open_dataset(base_dir / 'domain.lnd.CH_1km_navy.210407_new.nc')
surf_025deg = xr.open_dataset(base_dir / 'surfdata_CH_025deg_v2_hist_16pfts_Irrig_CMIP6_simyr2000_c230113.nc')
surf_05deg = xr.open_dataset(base_dir / 'surfdata_CH_05deg_v2_hist_16pfts_Irrig_CMIP6_simyr2000_c230113.nc')

input_global = xr.open_dataset(base_dir / 'surfdata_1km_CH_v3_hist_16pfts_Irrig_CMIP6_simyr2000_std25.nc')
input_highres = xr.open_dataset(base_dir / 'surfdata_1km_CH_v3_hist_16pfts_Irrig_CMIP6_NEW.nc')
lat_geo = input_global.LATIXY[:, 0]
lon_geo = input_global.LONGXY[0, :]
input_global = input_global.assign_coords({'lsmlon': lon_geo.data, 'lsmlat': lat_geo.data})
input_highres = input_highres.assign_coords({'lsmlon': lon_geo.data, 'lsmlat': lat_geo.data})
veg_global = input_global.PCT_NATVEG.where(domain.mask.data, np.nan)
veg_highres = input_highres.PCT_NATVEG.where(domain.mask.data, np.nan)

input_crujra = xr.open_dataset(base_dir / 'clmforc.crujra.TQ.2018-05.nc')
month_avg_cru = input_crujra.mean(dim='time')
del input_crujra
input_crujra_plus = xr.open_dataset(base_dir / 'clmforc.crujra.TQ.2018-05_lapse.nc')
month_avg_cruplus = input_crujra_plus.mean(dim='time')
month_avg_cruplus1 = month_avg_cruplus.TBOT.where(domain.mask.data, month_avg_cru.TBOT)
month_avg_cruplus['TBOT'] = month_avg_cruplus1
del input_crujra_plus
input_oshd = xr.open_dataset(base_dir / 'clmforc.OSHD1km.TQ.2018-05.nc')
month_avg_oshd = input_oshd.mean(dim='time')
month_avg_oshd1 = month_avg_oshd.TBOT.where(domain.mask.data, month_avg_cru.TBOT)
month_avg_oshd['TBOT'] = month_avg_oshd1
del input_oshd

lat_geo = month_avg_cruplus.LATIXY[:, 0]
lon_geo = month_avg_cruplus.LONGXY[0, :]

# create target grids for 0.25deg and 1deg
lat_new = surf_025deg.LATIXY.data[:, 0]
lon_new = surf_025deg.LONGXY.data[0, :]
ds_target_025 = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

lat_new = surf_05deg.LATIXY.data[:, 0]
lon_new = surf_05deg.LONGXY.data[0, :]
ds_target_05 = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

bf_out = base_dir / 'out_for_plots'
bf_out.mkdir(parents=True, exist_ok=True)

# output data, so I can work with it on my windows machine
file_oshd_05deg = bf_out / 'temp_oshd_05deg.nc'
file_cru_05deg = bf_out / 'temp_cru_05deg.nc'
file_cruplus_05deg = bf_out / 'temp_cruplus_05deg.nc'
file_oshd_025deg = bf_out / 'temp_oshd_025deg.nc'
file_cru_025deg = bf_out / 'temp_cru_025deg.nc'
file_cruplus_025deg = bf_out / 'temp_cruplus_025deg.nc'
file_oshd_1km = bf_out / 'temp_oshd_1km.nc'
file_cru_1km = bf_out / 'temp_cru_1km.nc'
file_cruplus_1km = bf_out / 'temp_cruplus_1km.nc'

file_all = [file_oshd_05deg, file_oshd_025deg, file_oshd_1km, file_cru_05deg, file_cru_025deg, file_cru_1km,
            file_cruplus_05deg, file_cruplus_025deg, file_cruplus_1km]

for file_check in file_all:
    if file_check.exists():
        print(file_check)
        os.remove(str(file_check))  #.unlink()

def regrid_data(ds_in, ds_target):
    ds_in1 = ds_in.assign_coords({'lsmlat': lat_geo, 'lsmlon': lon_geo})
    ds_in2 = ds_in1.swap_dims({'lon': 'lsmlon', 'lat': 'lsmlat'})
    ds_in3 = ds_in2.rename({'lsmlon': 'lon', 'lsmlat': 'lat'})

    regridder = xe.Regridder(ds_in3, ds_target, "nearest_s2d")  # nearest_s2d
    back_regridder = xe.Regridder(ds_target, ds_in3, "nearest_s2d")  # nearest_s2d
    ds_out1 = regridder(ds_in3)
    ds_out = ds_out1.where(ds_out1.TBOT != 0, np.nan)
    return ds_out, regridder, back_regridder


def regrid_oshd(ds_in, regridder, back_regridder):
    ds_in1 = ds_in.assign_coords({'lsmlat': lat_geo, 'lsmlon': lon_geo})
    ds_in2 = ds_in1.swap_dims({'lon': 'lsmlon', 'lat': 'lsmlat'})
    ds_in3 = ds_in2.rename({'lsmlon': 'lon', 'lsmlat': 'lat'})

    ds_out1 = regridder(ds_in3)
    ds_out = ds_out1.where(ds_out1.TBOT != 0, np.nan)

    month_avg_back1km = back_regridder(ds_out)
    ds_out_back = month_avg_back1km.TBOT.where(domain.mask.data, np.nan)

    return ds_out_back


month_avg_05deg, re, back_re = regrid_data(month_avg_cru, ds_target_05)
month_avg_05deg_back1km = back_re(month_avg_05deg)
temp_cru_05deg = month_avg_05deg_back1km.TBOT.where(domain.mask.data, np.nan)
temp_cru_05deg.to_netcdf(file_cru_05deg)

month_avg_05deg, re_05, back_re_05 = regrid_data(month_avg_cruplus, ds_target_05)
month_avg_05deg_back1km = back_re_05(month_avg_05deg)
temp_cruplus_05deg = month_avg_05deg_back1km.TBOT.where(domain.mask.data, np.nan)
temp_cruplus_05deg.to_netcdf(file_cruplus_05deg)

temp_oshd_05deg = regrid_oshd(month_avg_oshd, re_05, back_re_05)
temp_oshd_05deg.to_netcdf(file_oshd_05deg)

month_avg_025deg, re, back_re = regrid_data(month_avg_cru, ds_target_025)
month_avg_025deg_back1km = back_re(month_avg_025deg)
temp_cru_025deg = month_avg_025deg_back1km.TBOT.where(domain.mask.data, np.nan)
temp_cru_025deg.to_netcdf(file_cru_025deg)

month_avg_025deg, re_025, back_re_025 = regrid_data(month_avg_cruplus, ds_target_025)
month_avg_025deg_back1km = back_re_025(month_avg_025deg)
temp_cruplus_025deg = month_avg_025deg_back1km.TBOT.where(domain.mask.data, np.nan)
temp_cruplus_025deg.to_netcdf(file_cruplus_025deg)

temp_oshd_025deg = regrid_oshd(month_avg_oshd, re_025, back_re_025)
temp_oshd_025deg.to_netcdf(file_oshd_025deg)

temp_oshd_1km = month_avg_oshd.TBOT.where(domain.mask.data, np.nan)
temp_oshd_1km = temp_oshd_1km.assign_coords({'geolon': lon_geo, 'geolat': lat_geo})
ds1 = temp_oshd_1km.to_dataset()
ds2 = ds1.swap_dims({"lat": "geolat"})
temp_oshd_1km = ds2.swap_dims({"lon": "geolon"})
temp_oshd_1km.to_netcdf(file_oshd_1km)

temp_cru_1km = month_avg_cru.TBOT.where(domain.mask.data, np.nan)
temp_cru_1km = temp_cru_1km.assign_coords({'geolon': lon_geo, 'geolat': lat_geo})
ds1 = temp_cru_1km.to_dataset()
ds2 = ds1.swap_dims({"lat": "geolat"})
temp_cru_1km = ds2.swap_dims({"lon": "geolon"})
temp_cru_1km.to_netcdf(file_cru_1km)

temp_cruplus_1km = month_avg_cruplus.TBOT.where(domain.mask.data, np.nan)
temp_cruplus_1km = temp_cruplus_1km.assign_coords({'geolon': lon_geo, 'geolat': lat_geo})
ds1 = temp_cruplus_1km.to_dataset()
ds2 = ds1.swap_dims({"lat": "geolat"})
temp_cruplus_1km = ds2.swap_dims({"lon": "geolon"})
temp_cruplus_1km.to_netcdf(file_cruplus_1km)

# now vegetation
file_vegHR_025deg = bf_out / 'veg_hr_025.nc'
file_vegHR_05deg = bf_out / 'veg_hr_05.nc'
file_vegG_025deg = bf_out / 'veg_gl_025.nc'
file_vegG_05deg = bf_out / 'veg_gl_05.nc'
file_vegHR_1km = bf_out / 'veg_hr_1km.nc'
file_vegG_1km = bf_out / 'veg_gl_1km.nc'

file_all = [file_vegHR_025deg, file_vegHR_05deg, file_vegG_025deg, file_vegG_05deg, file_vegHR_1km, file_vegG_1km]

for file_check in file_all:
    if file_check.exists():
        file_check.unlink()

veg_global.to_netcdf(file_vegG_1km)
veg_highres.to_netcdf(file_vegHR_1km)

veg_highres_025deg_coarse = re_025(veg_highres)
veg_highres_05deg_coarse = re_05(veg_highres)
veg_global_025deg_coarse = re_025(input_global.PCT_NATVEG)
veg_global_05deg_coarse = re_05(input_global.PCT_NATVEG)

veg_highres_025deg_back = back_re_025(veg_highres_025deg_coarse)
veg_highres_05deg_back = back_re_05(veg_highres_05deg_coarse)
veg_highres_025deg = veg_highres_025deg_back.where(domain.mask.data, np.nan)
veg_highres_05deg = veg_highres_05deg_back.where(domain.mask.data, np.nan)

veg_global_025deg_back = back_re_025(veg_global_025deg_coarse)
veg_global_05deg_back = back_re_05(veg_global_05deg_coarse)
veg_global_025deg = veg_global_025deg_back.where(domain.mask.data, np.nan)
veg_global_05deg = veg_global_05deg_back.where(domain.mask.data, np.nan)

veg_highres_05deg.to_netcdf(file_vegHR_05deg)
veg_global_05deg.to_netcdf(file_vegG_05deg)
veg_global_025deg.to_netcdf(file_vegG_025deg)
veg_highres_025deg.to_netcdf(file_vegHR_025deg)
