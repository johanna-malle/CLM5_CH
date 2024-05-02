# -*- coding: utf-8 -*-
"""
Desc: for revisions we want to analyse everything at 0.25 and/or 0.5 degrees -> upscale 1km simulations
-> use xesmf_env_new2 !!!
Created on 16.01.24 10:33
@author: malle
"""

import numpy as np
from pathlib import Path
import xarray as xr
import rioxarray
import matplotlib
import cftime
import glob
import xesmf as xe
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cf


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# make switch to windows if working from home
mach = 'linux'
if mach == 'linux':
    bf = Path('/home/lud11/malle/CLM5_CH')
else:
    bf = Path('L:\malle\CLM5_CH')

# import surfdata & create target grid for xesmf:
surf_in = xr.open_dataset(bf / 'surfdata_1km_CH_v3_hist_16pfts_Irrig_CMIP6_NEW.nc')
lat_new = surf_in.LATIXY.data[:, 0]
lon_new = surf_in.LONGXY.data[0, :]
ds_target_1km = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

surf_in = xr.open_dataset(bf / 'surfdata_CH_025deg_v2_hist_16pfts_Irrig_CMIP6_simyr2000_c230113.nc')
lat_new = surf_in.LATIXY.data[:, 0]
lon_new = surf_in.LONGXY.data[0, :]
ds_target_025 = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

surf_in = xr.open_dataset(bf / 'surfdata_CH_05deg_v2_hist_16pfts_Irrig_CMIP6_simyr2000_c230113.nc')
lat_new = surf_in.LATIXY.data[:, 0]
lon_new = surf_in.LONGXY.data[0, :]
ds_target_05 = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

regridder_05 = xe.Regridder(ds_target_1km, ds_target_05, "conservative")
regridder_05_normed = xe.Regridder(ds_target_1km, ds_target_05, "conservative_normed")
regridder_05_bilinear = xe.Regridder(ds_target_1km, ds_target_05, "bilinear")
regridder_05_bilinear_extr_near = xe.Regridder(ds_target_1km, ds_target_05, "bilinear", extrap_method="nearest_s2d")
#regridder_05_bilinear_extr_ivd = xe.Regridder(ds_target_1km, ds_target_05, "bilinear", extrap_method="inverse_dist")


regridder_025 = xe.Regridder(ds_target_1km, ds_target_025, "conservative")
regridder_025_normed = xe.Regridder(ds_target_1km, ds_target_025, "conservative_normed")
regridder_025_bilinear = xe.Regridder(ds_target_1km, ds_target_025, "bilinear", extrap_method="nearest_s2d")
regridder_025_bilinear_extr_near = xe.Regridder(ds_target_1km, ds_target_025, "bilinear", extrap_method="nearest_s2d")
#regridder_025_bilinear_extr_ivd = xe.Regridder(ds_target_1km, ds_target_025, "bilinear", extrap_method="inverse_dist")


run_in_all = ['OSHD_FILES', 'OSHD_FILES_OLD', 'CRUJRA_FILES_noLapse', 'CRUJRA_FILES_noLapse_OLD', 'CRUJRA_FILES',
              'CRUJRA_FILES_OLD']

#  'OSHD_FILES_025_new', 'OSHD_FILES_05_new', 'CRUJRA_FILES_025deg_cru_new','CRUJRA_FILES_05deg_cru_new'
# 'CRUJRA_FILES_025deg_cru_new_lapse', 'CRUJRA_FILES_05deg_cru_new_lapse'

comp_025_OSHD = xr.open_dataset(glob.glob(str(bf / 'OSHD_FILES_025_new') + '/SNOW_DEPTH*')[0])
comp_05_OSHD = xr.open_dataset(glob.glob(str(bf / 'OSHD_FILES_05_new') + '/SNOW_DEPTH*')[0])

comp_025_CRU = xr.open_dataset(glob.glob(str(bf / 'CRUJRA_FILES_025deg_cru_new') + '/SNOW_DEPTH*')[0])
comp_05_CRU = xr.open_dataset(glob.glob(str(bf / 'CRUJRA_FILES_05deg_cru_new') + '/SNOW_DEPTH*')[0])

comp_025_CRU_LAPSE = xr.open_dataset(glob.glob(str(bf / 'CRUJRA_FILES_025deg_cru_new_lapse') + '/SNOW_DEPTH*')[0])
comp_05_CRU_LAPSE = xr.open_dataset(glob.glob(str(bf / 'CRUJRA_FILES_05deg_cru_new_lapse') + '/SNOW_DEPTH*')[0])

comp_all_05 = [comp_05_OSHD, comp_05_OSHD, comp_05_CRU, comp_05_CRU, comp_05_CRU_LAPSE, comp_05_CRU_LAPSE]
comp_all_025 = [comp_025_OSHD, comp_025_OSHD, comp_025_CRU, comp_025_CRU, comp_025_CRU_LAPSE, comp_025_CRU_LAPSE]

time_in = 1185
# now loop through all 1km simulations...
for id_run_in in range(np.size(run_in_all)):
    run_in = run_in_all[id_run_in]
    comp_05 = comp_all_05[id_run_in]
    comp_025 = comp_all_025[id_run_in]

    snow_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/SNOW_DEPTH*')[0])
    evap_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/QFLX_EVAP_TOT*')[0]) * 3600 * 24  # mm/s to mm/day
    evap_in_sum_yr = evap_in.resample(time='1Y').sum()

    snow_out_025 = regridder_025(snow_in)
    snow_out_025_normed = regridder_025_normed(snow_in)
    snow_out_025_bilinear = regridder_025_bilinear(snow_in)
    snow_out_025_bilinear_extr_near = regridder_025_bilinear_extr_near(snow_in)

    file_snow_025 = bf / run_in / 'upscale_025_SNOW_DEPTH.nc'
    file_snow_025_normed = bf / run_in / 'upscale_025_SNOW_DEPTH_normed.nc'
    file_snow_025_bilinear = bf / run_in / 'upscale_025_SNOW_DEPTH_bilinear.nc'
    file_snow_025_bilinear_ext = bf / run_in / 'upscale_025_SNOW_DEPTH_bilinear_ext.nc'

    snow_out_025.to_netcdf(file_snow_025)
    snow_out_025_normed.to_netcdf(file_snow_025_normed)
    snow_out_025_bilinear.to_netcdf(file_snow_025_bilinear)
    snow_out_025_bilinear_extr_near.to_netcdf(file_snow_025_bilinear_ext)

    snow_out_05 = regridder_05(snow_in)
    snow_out_05_normed = regridder_05_normed(snow_in)
    snow_out_05_bilinear = regridder_05_bilinear(snow_in)
    snow_out_05_bilinear_extr_near = regridder_05_bilinear_extr_near(snow_in)

    file_snow_05 = bf / run_in / 'upscale_05_SNOW_DEPTH.nc'
    file_snow_05_normed = bf / run_in / 'upscale_05_SNOW_DEPTH_normed.nc'
    file_snow_05_bilinear = bf / run_in / 'upscale_05_SNOW_DEPTH_bilinear.nc'
    file_snow_05_bilinear_ext = bf / run_in / 'upscale_05_SNOW_DEPTH_bilinear_ext.nc'

    snow_out_05.to_netcdf(file_snow_05)
    snow_out_05_normed.to_netcdf(file_snow_05_normed)
    snow_out_05_bilinear.to_netcdf(file_snow_05_bilinear)
    snow_out_05_bilinear_extr_near.to_netcdf(file_snow_05_bilinear_ext)

    font = {'size': 15}
    plt.rc('font', **font)
    proj_data = ccrs.PlateCarree()
    proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)

    a = np.array([[np.nanmax(snow_in.isel(time=time_in).DATA.data)],
                  [np.nanmax(snow_out_025.isel(time=time_in).DATA.data)],
                  [np.nanmax(snow_out_05.isel(time=time_in).DATA.data)]])
    min_all, max_all = np.nanmin(a), np.nanmax(a)

    # plot first now bilinear exp.
    fig, axs = plt.subplots(2, 3, frameon=False, figsize=[10, 6], subplot_kw={'projection': proj_map}, constrained_layout=True)
    plt.suptitle('Upscaling example for '+run_in+ ' : April 1st 2013')

    snow_in.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[0, 0], vmin=0, vmax=max_all)
    axs[0, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[0, 0].set_title('1km simulation')
    axs[0, 0].axis('off')
    axs[0, 0].axis('off')

    snow_out_025.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[0, 1], vmin=0, vmax=max_all)
    axs[0, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[0, 1].set_title('0.25 up, conserve')
    axs[0, 1].axis('off')
    axs[0, 1].axis('off')

    snow_out_05.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[0, 2], vmin=0, vmax=max_all)
    axs[0, 2].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[0, 2].set_title('0.5 up, conserve')
    axs[0, 2].axis('off')
    axs[0, 2].axis('off')

    comp_025.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[1, 1], vmin=0, vmax=max_all)
    axs[1, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[1, 1].set_title('0.25deg simulation')
    axs[1, 1].axis('off')
    axs[1, 1].axis('off')

    t1 = comp_05.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[1, 2], vmin=0, vmax=max_all)
    axs[1, 2].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[1, 2].set_title('0.5deg simulation')
    axs[1, 2].axis('off')
    axs[1, 2].axis('off')

    axs[1, 0].axis('off')
    axs[1, 0].axis('off')

    fig.colorbar(t1, ax=axs)
    plt.show()
    fig.savefig(bf / run_in / 'comp_conservative_regridding.png', facecolor='white', transparent=False)

    # plot first conservative comparison:
    fig, axs = plt.subplots(2, 3, frameon=False, figsize=[10, 6], subplot_kw={'projection': proj_map}, constrained_layout=True)
    plt.suptitle('Upscaling example for '+run_in+ ' : April 1st 2013')

    snow_in.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[0, 0], vmin=0, vmax=max_all)
    axs[0, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[0, 0].set_title('1km simulation')
    axs[0, 0].axis('off')
    axs[0, 0].axis('off')

    snow_out_025_bilinear.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[0, 1], vmin=0, vmax=max_all)
    axs[0, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[0, 1].set_title('0.25 up, bilinear')
    axs[0, 1].axis('off')
    axs[0, 1].axis('off')

    snow_out_05_bilinear.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[0, 2], vmin=0, vmax=max_all)
    axs[0, 2].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[0, 2].set_title('0.5 up, bilinear')
    axs[0, 2].axis('off')
    axs[0, 2].axis('off')

    comp_025.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[1, 1], vmin=0, vmax=max_all)
    axs[1, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[1, 1].set_title('0.25deg simulation')
    axs[1, 1].axis('off')
    axs[1, 1].axis('off')

    t1 = comp_05.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[1, 2], vmin=0, vmax=max_all)
    axs[1, 2].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[1, 2].set_title('0.5deg simulation')
    axs[1, 2].axis('off')
    axs[1, 2].axis('off')

    axs[1, 0].axis('off')
    axs[1, 0].axis('off')

    fig.colorbar(t1, ax=axs)
    plt.show()
    fig.savefig(bf / run_in / 'comp_bilinear_regridding.png', facecolor='white', transparent=False)

    # plot first now bilinear exp.
    fig, axs = plt.subplots(2, 3, frameon=False, figsize=[10, 6], subplot_kw={'projection': proj_map}, constrained_layout=True)
    plt.suptitle('Upscaling example for '+run_in+ ' : April 1st 2013')

    snow_in.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[0, 0], vmin=0, vmax=max_all)
    axs[0, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[0, 0].set_title('1km simulation')
    axs[0, 0].axis('off')
    axs[0, 0].axis('off')

    snow_out_025_bilinear_extr_near.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[0, 1], vmin=0, vmax=max_all)
    axs[0, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[0, 1].set_title('0.25 up, bilinear')
    axs[0, 1].axis('off')
    axs[0, 1].axis('off')

    snow_out_05_bilinear_extr_near.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[0, 2], vmin=0, vmax=max_all)
    axs[0, 2].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[0, 2].set_title('0.5 up, bilinear')
    axs[0, 2].axis('off')
    axs[0, 2].axis('off')

    comp_025.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[1, 1], vmin=0, vmax=max_all)
    axs[1, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[1, 1].set_title('0.25deg simulation')
    axs[1, 1].axis('off')
    axs[1, 1].axis('off')

    t1 = comp_05.isel(time=time_in).DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[1, 2], vmin=0, vmax=max_all)
    axs[1, 2].add_feature(cf.BORDERS, linewidth=1, edgecolor='darkred', alpha=1)
    axs[1, 2].set_title('0.5deg simulation')
    axs[1, 2].axis('off')
    axs[1, 2].axis('off')

    axs[1, 0].axis('off')
    axs[1, 0].axis('off')

    fig.colorbar(t1, ax=axs)
    plt.show()
    fig.savefig(bf / run_in / 'comp_bilinear_exp_regridding.png', facecolor='white', transparent=False)

    # now 0.25degree ET
    evap_out_025 = regridder_025(evap_in_sum_yr)
    evap_out_025_normed = regridder_025_normed(evap_in_sum_yr)
    evap_out_025_bilinear = regridder_025_bilinear(evap_in_sum_yr)
    evap_out_025_bilinear_ext = regridder_025_bilinear_extr_near(evap_in_sum_yr)

    file_025 = bf / run_in / 'upscale_025_QFLX_EVAP_TOT_YR_SUM.nc'
    file_025_normed = bf / run_in / 'upscale_025_QFLX_EVAP_TOT_YR_SUM_normed.nc'
    file_025_bilinear = bf / run_in / 'upscale_025_QFLX_EVAP_TOT_YR_SUM_bilinear.nc'
    file_025_bilinear_ext = bf / run_in / 'upscale_025_QFLX_EVAP_TOT_YR_SUM_bilinear_ext.nc'

    evap_out_025.to_netcdf(file_025)
    evap_out_025_normed.to_netcdf(file_025_normed)
    evap_out_025_bilinear.to_netcdf(file_025_bilinear)
    evap_out_025_bilinear_ext.to_netcdf(file_025_bilinear_ext)

    # now 0.5degree
    evap_out_05 = regridder_05(evap_in_sum_yr)
    evap_out_05_normed = regridder_05_normed(evap_in_sum_yr)
    evap_out_05_bilinear = regridder_05_bilinear(evap_in_sum_yr)
    evap_out_05_bilinear_ext = regridder_05_bilinear_extr_near(evap_in_sum_yr)

    file_05 = bf / run_in / 'upscale_05_QFLX_EVAP_TOT_YR_SUM.nc'
    file_05_normed = bf / run_in / 'upscale_05_QFLX_EVAP_TOT_YR_SUM_normed.nc'
    file_05_bilinear = bf / run_in / 'upscale_05_QFLX_EVAP_TOT_YR_SUM_bilinear.nc'
    file_05_bilinear_ext = bf / run_in / 'upscale_05_QFLX_EVAP_TOT_YR_SUM_bilinear_ext.nc'

    evap_out_05.to_netcdf(file_05)
    evap_out_05_normed.to_netcdf(file_05_normed)
    evap_out_05_bilinear.to_netcdf(file_05_bilinear)
    evap_out_05_bilinear_ext.to_netcdf(file_05_bilinear_ext)
