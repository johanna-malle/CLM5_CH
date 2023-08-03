# -*- coding: utf-8 -*-
"""
Desc:
Created on 25.01.23 15:23
@author: malle
"""

import numpy as np
from pathlib import Path
import xarray as xr
import sys

surfdata_new = xr.open_dataset('/media/malle/LaCie1/CLM5_input/surfdata_1km_CH_v3_hist_16pfts_Irrig_CMIP6_NEW.nc')
bf_out_new = Path('/home/malle/CLM5_install/input_data/PTCLM5_fluxnet_new')
bf_out_new_realPFT = Path('/home/malle/CLM5_install/input_data/PTCLM5_fluxnet_new_realPFT')

surfdata_orig = xr.open_dataset('/media/malle/LaCie1/CLM5_input/surfdata_1km_CH_v3_hist_16pfts_Irrig_CMIP6_NEW.nc')
bf_out_orig = Path('/home/malle/CLM5_install/input_data/PTCLM5_fluxnet_orig')
bf_out_orig_realPFT = Path('/home/malle/CLM5_install/input_data/PTCLM5_fluxnet_orig_realPFT')

name_locs = ['CH-Aws', 'CH-Cha', 'CH-Dav', 'CH-Fru', 'CH-Lae', 'CH-Oe2']
name_long = ['Alp Weissenstein', 'Chamau', 'Davos', 'Früebüel', 'Laegern', 'Oensingen crop']
lat_all = [46.5833, 47.2102, 46.8153, 47.1158, 47.4783, 47.2864]
lon_all = [9.790417, 8.4104, 9.8559, 8.5378, 8.3644, 7.7337]

for num_in in range(len(name_locs)):
    name_in = name_locs[num_in]
    name_long_in = name_long[num_in]
    lat_in = lat_all[num_in]
    lon_in = lon_all[num_in]

    file_out_new = bf_out_new / Path('surfdata_'+name_in+'.nc')
    file_out_orig = bf_out_orig / Path('surfdata_'+name_in+'.nc')
    file_out_new_realPFT = bf_out_new_realPFT / Path('surfdata_'+name_in+'.nc')
    file_out_orig_realPFT = bf_out_orig_realPFT / Path('surfdata_'+name_in+'.nc')

    close_lon = np.where(np.abs(surfdata_new.LONGXY[0, :]-lon_in) ==
                         np.min(np.abs(surfdata_new.LONGXY[0, :]-lon_in)))[0][0]
    close_lat = np.where(np.abs(surfdata_new.LATIXY[:, 0]-lat_in) ==
                         np.min(np.abs(surfdata_new.LATIXY[:, 0]-lat_in)))[0][0]

    surfdata_new_loc = surfdata_new.interp(lsmlon=[close_lon], lsmlat=[close_lat], method="nearest")
    surfdata_orig_loc = surfdata_orig.interp(lsmlon=[close_lon], lsmlat=[close_lat], method="nearest")
    surfdata_new_loc.LONGXY.values[:] = lon_in
    surfdata_new_loc.LATIXY.values[:] = lat_in

    surfdata_new_loc.to_netcdf(file_out_new)
    surfdata_orig_loc.to_netcdf(file_out_orig)

    surfdata_new_loc.LANDFRAC_PFT.values[:] = 1
    surfdata_new_loc.PCT_GLACIER.values[:] = 0
    surfdata_new_loc.PCT_LAKE.values[:] = 0
    surfdata_new_loc.PCT_WETLAND.values[:] = 0
    surfdata_new_loc.PCT_URBAN.values[:] = 0
    surfdata_new_loc.PCT_NATVEG.values[:] = 0
    surfdata_new_loc.PCT_CROP.values[:] = 0

    surfdata_orig_loc.LANDFRAC_PFT.values[:] = 1
    surfdata_orig_loc.PCT_GLACIER.values[:] = 0
    surfdata_orig_loc.PCT_LAKE.values[:] = 0
    surfdata_orig_loc.PCT_WETLAND.values[:] = 0
    surfdata_orig_loc.PCT_URBAN.values[:] = 0
    surfdata_orig_loc.PCT_NATVEG.values[:] = 0
    surfdata_orig_loc.PCT_CROP.values[:] = 0

    if name_in == 'CH-Oe2':
        surfdata_new_loc.PCT_CROP.values[:] = 100
        surfdata_new_loc.PCT_CFT.values[:] = 0
        surfdata_new_loc.PCT_CFT.values[0][:] = 100  # pft 15

        surfdata_orig_loc.PCT_CROP.values[:] = 100
        surfdata_orig_loc.PCT_CFT.values[:] = 0
        surfdata_orig_loc.PCT_CFT.values[0][:] = 100  # pft 15
    else:
        surfdata_new_loc.PCT_NAT_PFT.values[:] = 0
        surfdata_new_loc.PCT_NATVEG.values[:] = 100

        surfdata_orig_loc.PCT_NAT_PFT.values[:] = 0
        surfdata_orig_loc.PCT_NATVEG.values[:] = 100

    if name_in == 'CH-Aws':
        surfdata_new_loc.PCT_NAT_PFT.values[12][:] = 100  # 100% c3 arctic grass
        surfdata_orig_loc.PCT_NAT_PFT.values[12][:] = 100  # 100% c3 arctic grass
    elif name_in == 'CH-Cha':
        surfdata_new_loc.PCT_NAT_PFT.values[13][:] = 100  # 100% c3 grass
        surfdata_orig_loc.PCT_NAT_PFT.values[13][:] = 100  # 100% c3 grass
    elif name_in == 'CH-Dav':
        surfdata_new_loc.PCT_NAT_PFT.values[2][:] = 100  # 100% Needleleaf evergreen tree – boreal
        surfdata_orig_loc.PCT_NAT_PFT.values[2][:] = 100  # 100% Needleleaf evergreen tree – boreal
    elif name_in == 'CH-Fru':
        surfdata_new_loc.PCT_NAT_PFT.values[13][:] = 100  # 100% c3 gras
        surfdata_orig_loc.PCT_NAT_PFT.values[13][:] = 100  # 100% c3 gras
    elif name_in == 'CH-Lae':
        surfdata_new_loc.PCT_NAT_PFT.values[2][:] = 50  # 50% NET boreal
        surfdata_new_loc.PCT_NAT_PFT.values[8][:] = 50  # 50% BDT boreal

        surfdata_orig_loc.PCT_NAT_PFT.values[2][:] = 50  # 50% NET boreal
        surfdata_orig_loc.PCT_NAT_PFT.values[8][:] = 50  # 50% BDT boreal
    elif name_in == 'CH-Oe2':
        pass
    else:
        sys.exit()

    surfdata_new_loc.to_netcdf(file_out_new_realPFT)
    surfdata_orig_loc.to_netcdf(file_out_orig_realPFT)
