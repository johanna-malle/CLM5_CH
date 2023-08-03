# -*- coding: utf-8 -*-
"""
Desc: prep gpp, snow, ET files for further analysis/plots
Created on 13.02.23 10:54
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
id_surf_99 = surf_in.PCT_NATVEG > 99
id_surf_85 = surf_in.PCT_NATVEG > 85
id_surf_90 = surf_in.PCT_NATVEG > 90
lat_new = surf_in.LATIXY.data[:, 0]
lon_new = surf_in.LONGXY.data[0, :]
ds_target_1km = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

run_in_all = ['OSHD_FILES', 'OSHD_FILES_OLD', 'OSHD_FILES_025_new', 'OSHD_FILES_05_new',
              'CRUJRA_FILES_noLapse', 'CRUJRA_FILES_noLapse_OLD', 'CRUJRA_FILES_025deg_cru_new',
              'CRUJRA_FILES_05deg_cru_new', 'CRUJRA_FILES', 'CRUJRA_FILES_OLD',
              'CRUJRA_FILES_025deg_cru_new_lapse', 'CRUJRA_FILES_05deg_cru_new_lapse']

for run_in in run_in_all:
    gpp_in = xr.open_dataset(bf / run_in / 'FPSN_SP_MONTH_SUM.nc')  # already in units of gC m-2 mo-1
    snow_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/SNOW_DEPTH*')[0])
    fsno_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/FSNO*')[0])
    swe_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/H2OSNO*')[0])
    evap_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/QFLX_EVAP_TOT*')[0]) * 3600 * 24  # mm/s to mm/day
    evap_in_sum_yr = evap_in.resample(time='1Y').sum()
    if gpp_in.dims['lat'] != 272:
        file_int_gpp = bf / run_in / 're_FPSN_SP_MONTH_SUM.nc'
        file_int_snow = bf / run_in / 're_SNOW_DEPTH.nc'
        file_int_evap_yr = bf / run_in / 're_QFLX_EVAP_TOT_YR_SUM.nc'
        file_int_evap = bf / run_in / 're_QFLX_EVAP_TOT.nc'
        file_int_swe = bf / run_in / 're_H2OSNO.nc'
        file_int_fsno = bf / run_in / 're_FSNO.nc'
        if file_int_evap_yr.is_file():
            gpp_in = xr.open_dataset(file_int_gpp)
            snow_in = xr.open_dataset(file_int_snow)
            swe_in = xr.open_dataset(file_int_swe)
            evap_in_sum_yr = xr.open_dataset(file_int_evap_yr)
            evap_in = xr.open_dataset(file_int_evap)
            fsno_in = xr.open_dataset(file_int_fsno)
        else:
            gpp_in_coarse = xr.open_dataset(bf / run_in / 'FPSN_SP_MONTH_SUM.nc')
            regridder_nn = xe.Regridder(gpp_in_coarse, ds_target_1km, "nearest_s2d")
            gpp_in = regridder_nn(gpp_in_coarse)
            gpp_in.to_netcdf(file_int_gpp)

            snow_in_coarse = xr.open_dataset(bf / run_in / 'SNOW_DEPTH.nc')
            snow_in = regridder_nn(snow_in_coarse)
            snow_in.to_netcdf(file_int_snow)

            swe_in = regridder_nn(swe_in)
            swe_in.to_netcdf(file_int_swe)

            fsno_in = regridder_nn(fsno_in)
            fsno_in.to_netcdf(file_int_fsno)

            evap_in_sum_yr = regridder_nn(evap_in_sum_yr)
            evap_in_sum_yr.to_netcdf(file_int_evap_yr)

            evap_in = regridder_nn(evap_in)
            evap_in.to_netcdf(file_int_evap)

    datetimes1 = snow_in.indexes['time']
    datetimes1_evapo = evap_in_sum_yr.indexes['time']

    yr_all = [2016, 2017, 2018, 2019]
    print(run_in)
    for yr_in in yr_all:
        gpp_07 = gpp_in.DATA.sel(time=cftime._cftime.DatetimeNoLeap(yr_in, 7, 1), method='nearest')
        gpp_08 = gpp_in.DATA.sel(time=cftime._cftime.DatetimeNoLeap(yr_in, 8, 1), method='nearest')
        gpp_avg = (gpp_07 + gpp_08) / 2

        gpp_05 = gpp_in.DATA.sel(time=cftime._cftime.DatetimeNoLeap(yr_in, 5, 1), method='nearest')
        gpp_06 = gpp_in.DATA.sel(time=cftime._cftime.DatetimeNoLeap(yr_in, 6, 1), method='nearest')
        gpp_avg_0506 = (gpp_05 + gpp_06) / 2

        id_start = np.where(datetimes1 == cftime.DatetimeNoLeap(yr_in, 1, 1, 0, 0, 0, 0, has_year_zero=True))[0][0]
        id_end = np.where(datetimes1 == cftime.DatetimeNoLeap(yr_in, 6, 30, 0, 0, 0, 0, has_year_zero=True))[0][0]
        id_end_dec = np.where(datetimes1 == cftime.DatetimeNoLeap(yr_in, 12, 31, 0, 0, 0, 0, has_year_zero=True))[0][0]
        id_start_hydro = np.where(datetimes1 == cftime.DatetimeNoLeap(yr_in-1, 10, 1, 0, 0, 0, 0,
                                                                      has_year_zero=True))[0][0]
        id_end_hydro = np.where(datetimes1 == cftime.DatetimeNoLeap(yr_in, 9, 30, 0, 0, 0, 0, has_year_zero=True))[0][0]
        id_evapo = np.where(datetimes1_evapo.year == yr_in)
        evapo_yr = evap_in_sum_yr.isel(time=id_evapo[0][0])
        evapo_hydro = evap_in.loc[{'time': slice(datetimes1[id_start_hydro], datetimes1[id_end_hydro])}]
        evapo_hydro_sum = evapo_hydro.DATA.sum(dim='time')

        snow_jan_jun = snow_in.loc[{'time': slice(datetimes1[id_start], datetimes1[id_end])}]
        snow_jan_dec = snow_in.loc[{'time': slice(datetimes1[id_start], datetimes1[id_end_dec])}]
        snow_hydro = snow_in.loc[{'time': slice(datetimes1[id_start_hydro], datetimes1[id_end_hydro])}]
        swe_hydro = swe_in.loc[{'time': slice(datetimes1[id_start_hydro], datetimes1[id_end_hydro])}]

        # count number of snow days with more than 2cm snow on ground
        snow_count = snow_jan_jun.DATA.where(snow_jan_jun.DATA > 0.02).count(dim='time')
        snow_count_year = snow_jan_dec.DATA.where(snow_jan_dec.DATA > 0.02).count(dim='time')
        snow_count_hydroyear = snow_hydro.DATA.where(snow_hydro.DATA > 0.02).count(dim='time')

        # calculate total amount of snow during (a) jan-jun, (b) jan-dec, and (c) hydrological year
        snow_total = snow_jan_jun.DATA.sum(dim='time')
        snow_total_year = snow_jan_dec.DATA.sum(dim='time')
        snow_total_hydroyear = snow_hydro.DATA.sum(dim='time')
        swe_total_hydroyear = swe_hydro.DATA.sum(dim='time')
        swe_inc = swe_hydro.DATA.diff(dim='time')
        swe_inc = swe_inc.where(swe_inc > 0, np.nan)
        swe_inc_total = swe_inc.sum(dim='time')

        # calculate melt-out date: 1) find max. SWE, then first time after that where SWE == 0
        swe_hydro_nonan = swe_hydro.fillna(20)  # dummy way of dealing with NaNs
        id_peak_snow = swe_hydro_nonan.DATA.argmax(dim='time')
        cft_peak_snow = swe_hydro.time[id_peak_snow]
        datetime_peak_snow = np.datetime64(datetimes1[id_start_hydro].isoformat()) \
            + id_peak_snow.data.astype('timedelta64[D]')
        swe_hydro_post = swe_hydro.where(swe_hydro.time > cft_peak_snow, np.nan)

        melt_out_date = swe_hydro_post.DATA.where(swe_hydro_post['DATA'] == 0).\
            idxmin(dim='time').fillna(np.datetime64('now')).data.astype(str).astype(np.datetime64)

        melt_out_date_hydroyear = snow_count_hydroyear.copy()
        melt_out_date_hydroyear.data = melt_out_date

        melt_out_doy = np.zeros(melt_out_date.shape)
        for i in range(1, melt_out_date.shape[0]):
            dates_in = pd.to_datetime(melt_out_date[i, :]).to_series()
            doy_out = dates_in.dt.dayofyear
            melt_out_doy[i, :] = doy_out

        melt_out_doy_hydroyear = snow_count_hydroyear.copy()
        melt_out_doy_hydroyear.data = melt_out_doy

        # add masks for pixels<1000, 1000-2000, 2000-3000 as well as vegetation percentages for later analysis
        dem = rioxarray.open_rasterio(bf / 'BAFU_DEM_2020_1000.tif')
        id_1000 = (dem > 0) & (dem < 1000)
        id_2000 = (dem >= 1000) & (dem < 2000)
        id_3000 = (dem >= 2000)

        snow_count.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        snow_count_year.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        snow_count_hydroyear.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        gpp_avg.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        gpp_avg_0506.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        evapo_yr.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        evapo_hydro_sum.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        snow_total.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        snow_total_year.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        snow_total_hydroyear.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        melt_out_date_hydroyear.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        melt_out_doy_hydroyear.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        swe_total_hydroyear.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)
        swe_inc_total.coords['mask_veg'] = (('lat', 'lon'), (id_surf_99 * 1).squeeze().data)

        snow_count.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        snow_count_year.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        snow_count_hydroyear.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        gpp_avg.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        gpp_avg_0506.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        evapo_yr.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        evapo_hydro_sum.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        snow_total.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        snow_total_year.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        snow_total_hydroyear.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        melt_out_date_hydroyear.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        melt_out_doy_hydroyear.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        swe_total_hydroyear.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)
        swe_inc_total.coords['mask_veg_85'] = (('lat', 'lon'), (id_surf_85 * 1).squeeze().data)

        snow_count.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        snow_count_year.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        snow_count_hydroyear.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        gpp_avg.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        gpp_avg_0506.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        evapo_yr.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        evapo_hydro_sum.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        snow_total.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        snow_total_year.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        snow_total_hydroyear.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        melt_out_date_hydroyear.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        melt_out_doy_hydroyear.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        swe_total_hydroyear.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)
        swe_inc_total.coords['mask_veg_90'] = (('lat', 'lon'), (id_surf_90 * 1).squeeze().data)

        snow_count.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        snow_count_year.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        snow_count_hydroyear.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        gpp_avg.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        gpp_avg_0506.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        evapo_yr.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        evapo_hydro_sum.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        snow_total.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        snow_total_year.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        snow_total_hydroyear.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        melt_out_date_hydroyear.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        melt_out_doy_hydroyear.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        swe_total_hydroyear.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))
        swe_inc_total.coords['mask_1000'] = (('lat', 'lon'), np.flipud((id_1000*1).squeeze().data))

        snow_count.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        snow_count_year.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        snow_count_hydroyear.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        gpp_avg.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        gpp_avg_0506.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        evapo_yr.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        evapo_hydro_sum.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        snow_total.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        snow_total_year.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        snow_total_hydroyear.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        melt_out_date_hydroyear.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        melt_out_doy_hydroyear.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        swe_total_hydroyear.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))
        swe_inc_total.coords['mask_2000'] = (('lat', 'lon'), np.flipud((id_2000*1).squeeze().data))

        snow_count.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        snow_count_year.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        snow_count_hydroyear.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        gpp_avg.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        gpp_avg_0506.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        evapo_yr.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        evapo_hydro_sum.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        snow_total.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        snow_total_year.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        snow_total_hydroyear.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        melt_out_date_hydroyear.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        melt_out_doy_hydroyear.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        swe_total_hydroyear.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))
        swe_inc_total.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_3000*1).squeeze().data))

        # set to nan what is outside of domaain
        snow_count1 = snow_count.where(~np.isnan(gpp_avg), np.nan)
        snow_count_yr1 = snow_count_year.where(~np.isnan(gpp_avg), np.nan)
        snow_count_hydroyr1 = snow_count_hydroyear.where(~np.isnan(gpp_avg), np.nan)
        evapo_yr1 = evapo_yr.where(~np.isnan(gpp_avg), np.nan)
        evapo_hydro_sum1 = evapo_hydro_sum.where(~np.isnan(gpp_avg), np.nan)
        snow_total1 = snow_total.where(~np.isnan(gpp_avg), np.nan)
        snow_total_yr1 = snow_total_year.where(~np.isnan(gpp_avg), np.nan)
        snow_total_hydroyear1 = snow_total_hydroyear.where(~np.isnan(gpp_avg), np.nan)
        melt_out_date_hydroyear1 = melt_out_date_hydroyear.where(~np.isnan(gpp_avg), np.datetime64('NaT'))
        melt_out_doy_hydroyear1 = melt_out_doy_hydroyear.where(~np.isnan(gpp_avg), np.nan)
        swe_total_hydroyear1 = swe_total_hydroyear.where(~np.isnan(gpp_avg), np.nan)
        swe_inc_total1 = swe_inc_total.where(~np.isnan(gpp_avg), np.nan)

        # write to netcdf files for later analysis
        ending = '_v0.nc'
        file_out_gpp = 'gpp_avg_' + str(yr_in) + ending
        file_out_gpp_0506 = 'gpp_avg_0506_' + str(yr_in) + ending
        file_out_snow = 'snow_count_' + str(yr_in) + ending
        file_out_snow_yr = 'snow_count_year_' + str(yr_in) + ending
        file_out_snow_hydroyr = 'snow_count_HYDROyear_' + str(yr_in) + ending
        file_out_evap_yr = 'total_ET_year_' + str(yr_in) + ending
        file_out_evap_hydroyr = 'total_ET_HYDORyear_' + str(yr_in) + ending
        file_out_snow_total = 'snow_total_' + str(yr_in) + ending
        file_out_snow_total_yr = 'snow_total_year_' + str(yr_in) + ending
        file_out_snow_total_hydroyr = 'snow_total_HYDROyear_' + str(yr_in) + ending
        file_out_melt_date = 'melt_out_' + str(yr_in) + ending
        file_out_melt_doy = 'melt_doy_' + str(yr_in) + ending
        file_out_swe_total = 'swe_total_HYDROyear_' + str(yr_in) + ending
        file_out_swe_inc_total = 'swe_inc_total_HYDROyear_' + str(yr_in) + ending

        gpp_avg.to_netcdf(bf / run_in / file_out_gpp)
        gpp_avg_0506.to_netcdf(bf / run_in / file_out_gpp_0506)
        snow_count1.to_netcdf(bf / run_in / file_out_snow)
        snow_count_yr1.to_netcdf(bf / run_in / file_out_snow_yr)
        snow_count_hydroyr1.to_netcdf(bf / run_in / file_out_snow_hydroyr)
        evapo_yr1.to_netcdf(bf / run_in / file_out_evap_yr)
        evapo_hydro_sum1.to_netcdf(bf / run_in / file_out_evap_hydroyr)
        snow_total1.to_netcdf(bf / run_in / file_out_snow_total)
        snow_total_yr1.to_netcdf(bf / run_in / file_out_snow_total_yr)
        snow_total_hydroyear1.to_netcdf(bf / run_in / file_out_snow_total_hydroyr)
        melt_out_date_hydroyear1.to_netcdf(bf / run_in / file_out_melt_date)
        melt_out_doy_hydroyear1.to_netcdf(bf / run_in / file_out_melt_doy)
        swe_total_hydroyear1.to_netcdf(bf / run_in / file_out_swe_total)
        swe_inc_total1.to_netcdf(bf / run_in / file_out_swe_inc_total)
