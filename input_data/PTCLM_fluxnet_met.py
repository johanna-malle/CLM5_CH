# -*- coding: utf-8 -*-
"""
Desc:
Created on 25.01.23 11:40
@author: malle
"""

import numpy as np
import pandas as pd
from pathlib import Path
import xarray as xr
import xesmf as xe
import glob

# first meteo data
bf_oshd = Path('/media/malle/LaCie1/CLM5_input/OSHD_CH1km_v3')
bf_cru = Path('/media/malle/LaCie1/CLM5_input/cru_jra_CLM5_CH')
bf_out = Path('/media/malle/LaCie1/CLM5_input/PTCLM/met_input_fluxnet')

vars_all = ['lwr_pa', 'precip', 'solar', 'tair_rh', 'wind_zbot']
time_all = pd.date_range(start='2013-01-01', end='2020-12-31', freq='M')
time_all_str = time_all.strftime('%Y-%m, %r')

name_locs = ['CH-Aws', 'CH-Cha', 'CH-Dav', 'CH-Fru', 'CH-Lae', 'CH-Oe2']
name_long = ['Alp Weissenstein', 'Chamau', 'Davos', 'Früebüel', 'Laegern', 'Oensingen crop']
lat_all = [46.5833, 47.2102, 46.8153, 47.1158, 47.4783, 47.2864]
lon_all = [9.790417, 8.4104, 9.8559, 8.5378, 8.3644, 7.7337]

# also need to regrid 1km to 0.5 and 0.25deg
bf = Path('/home/lud11/malle/CLM5_CH')
surf_025deg = xr.open_dataset(bf / 'surfdata_CH_025deg_v2_hist_16pfts_Irrig_CMIP6_simyr2000_c230113.nc')
surf_05deg = xr.open_dataset(bf / 'surfdata_CH_05deg_v2_hist_16pfts_Irrig_CMIP6_simyr2000_c230113.nc')

lat_new = surf_025deg.LATIXY.data[:, 0]
lon_new = surf_025deg.LONGXY.data[0, :]
ds_target_025 = xr.Dataset({"lat": (["lat"], lat_new),
                            "lon": (["lon"], lon_new)})
lat_new = surf_05deg.LATIXY.data[:, 0]
lon_new = surf_05deg.LONGXY.data[0, :]
ds_target_05 = xr.Dataset({"lat": (["lat"], lat_new),
                          "lon": (["lon"], lon_new)})

file_ref = xr.open_dataset('/media/malle/LaCie1/CLM5_input/surfdata_1km_CH_v3_hist_16pfts_Irrig_CMIP6_NEW.nc')
file_ref1 = file_ref.assign_coords({'lsmlat': file_ref.LATIXY[:, 0], 'lsmlon': file_ref.LONGXY[0, :]})
file_ref2 = file_ref1.rename({'lsmlon': 'lon', 'lsmlat': 'lat'})

regridder_nearest_05 = xe.Regridder(file_ref2, ds_target_05, "nearest_s2d")
regridder_nearest_025 = xe.Regridder(file_ref2, ds_target_025, "nearest_s2d")

for time_in in range(len(time_all_str)):
    date_in = time_all_str[time_in][:7]

    for var in vars_all:
        file_in_oshd = glob.glob(str(bf_oshd)+'/'+var+'/*'+date_in+'*')
        oshd_in = xr.open_dataset(file_in_oshd[0])
        oshd_in_05 = regridder_nearest_05(oshd_in)
        oshd_in_025 = regridder_nearest_025(oshd_in)

        file_in_cru = glob.glob(str(bf_cru)+'/'+var+'/*'+date_in+'*')
        cru_in = xr.open_dataset(file_in_cru[0])
        cru_in_05 = regridder_nearest_05(cru_in)
        cru_in_025 = regridder_nearest_025(cru_in)
        if var == 'tair_rh':
            file_in_cru_lapse = glob.glob(str(bf_cru) + '/' + var + '_lapse/*' + date_in + '*')
            cru_in_lapse = xr.open_dataset(file_in_cru_lapse[0])
            cru_in_05_lapse = regridder_nearest_05(cru_in_lapse)
            cru_in_025_lapse = regridder_nearest_025(cru_in_lapse)

        for num_in in range(len(name_locs)):
            name_in = name_locs[num_in]
            name_long_in = name_long[num_in]
            lat_in = lat_all[num_in]
            lon_in = lon_all[num_in]

            close_lon = np.where(np.abs(oshd_in.LONGXY[0, :] - (lon_in + 180)) == np.min(np.abs(oshd_in.LONGXY[0, :] - (lon_in + 180))))[0][0]
            close_lat = np.where(np.abs(oshd_in.LATIXY[:, 0] - lat_in) == np.min(np.abs(oshd_in.LATIXY[:, 0] - lat_in)))[0][0]
            oshd_loc = oshd_in.interp(lon=[close_lon], lat=[close_lat], method="nearest")

            oshd_loc_025 = oshd_in_025.interp(lon=[lon_in], lat=[lat_in], method="nearest")
            oshd_loc_05 = oshd_in_05.interp(lon=[lon_in], lat=[lat_in], method="nearest")

            oshd_out = bf_out / name_in / 'oshd' / var
            oshd_out.mkdir(parents=True, exist_ok=True)
            oshd_out_025 = bf_out / name_in / 'oshd_025' / var
            oshd_out_025.mkdir(parents=True, exist_ok=True)
            oshd_out_05 = bf_out / name_in / 'oshd_05' / var
            oshd_out_05.mkdir(parents=True, exist_ok=True)

            oshd_loc.to_netcdf(oshd_out / Path(file_in_oshd[0].split('/')[-1]))
            oshd_loc_025.to_netcdf(oshd_out_025 / Path(file_in_oshd[0].split('/')[-1]))
            oshd_loc_05.to_netcdf(oshd_out_05 / Path(file_in_oshd[0].split('/')[-1]))

            cru_loc = cru_in.interp(lon=[close_lon], lat=[close_lat], method="nearest")
            cru_loc_025 = cru_in_025.interp(lon=[lon_in], lat=[lat_in], method="nearest")
            cru_loc_05 = cru_in_05.interp(lon=[lon_in], lat=[lat_in], method="nearest")

            cru_out = bf_out / name_in / 'cru' / var
            cru_out.mkdir(parents=True, exist_ok=True)
            cru_out_025 = bf_out / name_in / 'cru_025' / var
            cru_out_025.mkdir(parents=True, exist_ok=True)
            cru_out_05 = bf_out / name_in / 'cru_05' / var
            cru_out_05.mkdir(parents=True, exist_ok=True)

            cru_loc.to_netcdf(cru_out / Path(file_in_cru[0].split('/')[-1]))
            cru_loc_025.to_netcdf(cru_out_025 / Path(file_in_cru[0].split('/')[-1]))
            cru_loc_05.to_netcdf(cru_out_05 / Path(file_in_cru[0].split('/')[-1]))

            if var == 'tair_rh':
                cru_lapse_loc = cru_in_lapse.interp(lon=[close_lon], lat=[close_lat], method="nearest")
                cru_lapse_loc_025 = cru_in_025_lapse.interp(lon=[lon_in], lat=[lat_in], method="nearest")
                cru_lapse_loc_05 = cru_in_05_lapse.interp(lon=[lon_in], lat=[lat_in], method="nearest")

                cru_out = bf_out / name_in / 'cru' / 'tair_rh_lapse'
                cru_out.mkdir(parents=True, exist_ok=True)
                cru_out_025 = bf_out / name_in / 'cru_025' / 'tair_rh_lapse'
                cru_out_025.mkdir(parents=True, exist_ok=True)
                cru_out_05 = bf_out / name_in / 'cru_05' / 'tair_rh_lapse'
                cru_out_05.mkdir(parents=True, exist_ok=True)

                cru_lapse_loc.to_netcdf(cru_out / Path(file_in_cru_lapse[0].split('/')[-1]))
                cru_lapse_loc_025.to_netcdf(cru_out_025 / Path(file_in_cru_lapse[0].split('/')[-1]))
                cru_lapse_loc_05.to_netcdf(cru_out_05 / Path(file_in_cru_lapse[0].split('/')[-1]))