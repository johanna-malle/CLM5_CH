# -*- coding: utf-8 -*-
"""
Desc:
Created on 08.12.22 18:02
@author: malle
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib
import cftime
import xesmf as xe
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import platform

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# make switch to windows if working from home
if platform.system() == 'Windows':
    bf = Path('L:\malle\CLM5_CH')
else:
    bf = Path('/home/lud11/malle/CLM5_CH')

g_carbon = 12.0107
conv_umol = g_carbon * 3600 * 24 * 0.000001  # fluxnet data is measured in umol/s -> convert...

surf_025deg = xr.open_dataset(bf / 'surfdata_CH_025deg_v2_hist_16pfts_Irrig_CMIP6_simyr2000_c230113.nc')
surf_05deg = xr.open_dataset(bf / 'surfdata_CH_05deg_v2_hist_16pfts_Irrig_CMIP6_simyr2000_c230113.nc')

lat_new = surf_025deg.LATIXY.data[:, 0]
lon_new = surf_025deg.LONGXY.data[0, :]
ds_target_025 = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

lat_new = surf_05deg.LATIXY.data[:, 0]
lon_new = surf_05deg.LONGXY.data[0, :]
ds_target_05 = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

# load clm5 simulations for the various model configurations -> focus on july/august average 2017
gpp_load = xr.open_dataset(bf / 'OSHD_FILES' / 'FPSN_SP_MONTH_SUM.nc')
gpp_201707 = gpp_load.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 7, 1), method='nearest') * conv_umol
gpp_201708 = gpp_load.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 8, 1), method='nearest') * conv_umol
gpp_oshd_avg = (gpp_201707 + gpp_201708) / 2

gpp_load = xr.open_dataset(bf / 'OSHD_FILES_OLD' / 'FPSN_SP_MONTH_SUM.nc')
gpp_201707 = gpp_load.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 7, 1), method='nearest') * conv_umol
gpp_201708 = gpp_load.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 8, 1), method='nearest') * conv_umol
gpp_oshd_global_avg = (gpp_201707 + gpp_201708) / 2

gpp_load = xr.open_dataset(bf / 'CRUJRA_FILES_OLD' / 'FPSN_SP_MONTH_SUM.nc')
gpp_201707 = gpp_load.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 7, 1), method='nearest') * conv_umol
gpp_201708 = gpp_load.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 8, 1), method='nearest') * conv_umol
gpp_crujra_avg = (gpp_201707 + gpp_201708) / 2

gpp_load = xr.open_dataset(bf / 'CRUJRA_FILES_noLapse_OLD' / 'FPSN_SP_MONTH_SUM.nc')
gpp_201707 = gpp_load.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 7, 1), method='nearest') * conv_umol
gpp_201708 = gpp_load.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 8, 1), method='nearest') * conv_umol
gpp_crujra_nolapse_avg = (gpp_201707 + gpp_201708) / 2

gpp_load = xr.open_dataset(bf / 'CRUJRA_FILES_025deg_cru_new' / 'FPSN.nc')  # not monthly sums yet!
sum_month = gpp_load.resample(time='1M').sum()
bacK_regridder_nearest_025 = xe.Regridder(ds_target_025, gpp_crujra_nolapse_avg, "nearest_s2d")
sum_month_re = bacK_regridder_nearest_025(sum_month)
gpp_201707 = sum_month_re.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 7, 1), method='nearest') * conv_umol
gpp_201708 = sum_month_re.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 8, 1), method='nearest') * conv_umol
gpp_crujra_025_avg = (gpp_201707 + gpp_201708) / 2

gpp_load = xr.open_dataset(bf / 'CRUJRA_FILES_05deg_cru_new' / 'FPSN.nc')
sum_month = gpp_load.resample(time='1M').sum()
bacK_regridder_nearest_05 = xe.Regridder(ds_target_05, gpp_crujra_nolapse_avg, "nearest_s2d")
sum_month_re = bacK_regridder_nearest_05(sum_month)
gpp_201707 = sum_month_re.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 7, 1), method='nearest') * conv_umol
gpp_201708 = sum_month_re.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 8, 1), method='nearest') * conv_umol
gpp_crujra_05_avg = (gpp_201707 + gpp_201708) / 2

# now coarse OSHD
gpp_load = xr.open_dataset(bf / 'OSHD_FILES_05_new' / 'FPSN.nc')
sum_month = gpp_load.resample(time='1M').sum()
sum_month_re = bacK_regridder_nearest_05(sum_month)
gpp_201707 = sum_month_re.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 7, 1), method='nearest') * conv_umol
gpp_201708 = sum_month_re.DATA.sel(time=cftime._cftime.DatetimeNoLeap(2017, 8, 1), method='nearest') * conv_umol
gpp_oshd_05_avg = (gpp_201707 + gpp_201708) / 2

surf_in = xr.open_dataset(bf / 'surfdata_1km_CH_v3_hist_16pfts_Irrig_CMIP6_NEW.nc')
id_surf = surf_in.PCT_LAKE < 101

# in case we want to mask out lakes:
gpp_crujra_025_avg.coords['mask_lake'] = (('lat', 'lon'), (id_surf*1).squeeze().data)
gpp_crujra_nolapse_avg.coords['mask_lake'] = (('lat', 'lon'), (id_surf*1).squeeze().data)
gpp_crujra_avg.coords['mask_lake'] = (('lat', 'lon'), (id_surf*1).squeeze().data)
gpp_oshd_global_avg.coords['mask_lake'] = (('lat', 'lon'), (id_surf*1).squeeze().data)
gpp_oshd_avg.coords['mask_lake'] = (('lat', 'lon'), (id_surf*1).squeeze().data)

id_mask_nan = gpp_crujra_025_avg != 0
gpp_crujra_025_avg.coords['mask_nan'] = (('lat', 'lon'), (id_mask_nan*1).squeeze().data)

id_mask_1km = ~np.isnan(gpp_oshd_avg)
gpp_crujra_025_avg.coords['mask_nan_1km'] = (('lat', 'lon'), (id_mask_1km*1).squeeze().data)

a = np.array([[gpp_oshd_avg.data], [gpp_oshd_global_avg.data], [gpp_crujra_nolapse_avg.data],
              [gpp_oshd_05_avg.data]])
a1 = np.array([[gpp_oshd_avg.data - gpp_oshd_global_avg.data], [gpp_crujra_nolapse_avg.data - gpp_oshd_global_avg.data],
              [gpp_oshd_05_avg.data - gpp_oshd_global_avg.data]])
min_gpp_all, max_gpp_all = np.nanmin(a), np.nanmax(a)
min_gpp_diff, max_gpp_diff = np.nanmin(a1), np.nanmax(a1)
diff_lim = np.max([[np.abs(min_gpp_diff)], [np.abs(max_gpp_diff)]])

# make plot:
font = {'size': 15}
plt.rc('font', **font)
proj_data = ccrs.PlateCarree()
proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)

fig, axs = plt.subplots(4, 2, figsize=[9, 15], frameon=False, subplot_kw={'projection': proj_map})
gpp_oshd_global_avg.where((gpp_oshd_avg.mask_lake == 1)).plot(transform=proj_data, cmap='viridis', add_colorbar=False,
                                                              vmin=min_gpp_all, vmax=max_gpp_all, ax=axs[0, 1])
axs[0, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
plt.title('')
axs[0, 1].text(-0.07, 0.55, 'Clim$_{OSHD 1km}$+LU$_{Gl 1km}$', va='bottom', ha='center', rotation='vertical',
               rotation_mode='anchor', transform=axs[0, 1].transAxes, fontsize=17)
axs[0, 1].axis('off')
axs[0, 0].axis('off')
# effect of land use info:
gpp_oshd_avg.where((gpp_oshd_avg.mask_lake == 1)).plot(transform=proj_data, cmap='viridis',
                                                       add_colorbar=False, vmin=min_gpp_all,
                                                       vmax=max_gpp_all, ax=axs[1, 0])
axs[1, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
plt.title('')
axs[1, 0].axis('off')
axs[1, 0].text(-0.07, 0.55, 'Clim$_{OSHD 1km}$+LU$_{HR 1km}$', va='bottom', ha='center', rotation='vertical',
               rotation_mode='anchor', transform=axs[1, 0].transAxes, fontsize=17)
(gpp_oshd_avg.where((gpp_oshd_avg.mask_lake == 1)) - gpp_oshd_global_avg.where((gpp_oshd_avg.mask_lake == 1))).\
    plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False, vmin=-diff_lim, vmax=diff_lim, ax=axs[1, 1])
axs[1, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
plt.title('')
axs[1, 1].axis('off')
# effect of climate forcing:
gpp_crujra_nolapse_avg.where((gpp_oshd_avg.mask_lake == 1)).plot(transform=proj_data, cmap='viridis',
                                                                 add_colorbar=False, vmin=min_gpp_all,
                                                                 vmax=max_gpp_all, ax=axs[2, 0])
axs[2, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
plt.title('')
axs[2, 0].axis('off')
axs[2, 0].text(-0.07, 0.55, 'Clim$_{CRU 1km}$+LU$_{Gl 1km}$', va='bottom', ha='center',
               rotation='vertical', rotation_mode='anchor', transform=axs[2, 0].transAxes, fontsize=17)
(gpp_crujra_nolapse_avg.where((gpp_oshd_avg.mask_lake == 1)) -
 gpp_oshd_global_avg.where((gpp_oshd_avg.mask_lake == 1))).plot(transform=proj_data, cmap='RdBu_r',
                                                                add_colorbar=False, vmin=-diff_lim,
                                                                vmax=diff_lim, ax=axs[2, 1])
axs[2, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
plt.title('')
axs[2, 1].axis('off')
# effect of resolution:
p = gpp_oshd_05_avg.where((gpp_oshd_avg.mask_lake == 1) & (gpp_crujra_025_avg.mask_nan == 1) &
                          (gpp_crujra_025_avg.mask_nan_1km == 1)).plot(transform=proj_data, cmap='viridis',
                                                                       add_colorbar=False, vmin=min_gpp_all,
                                                                       vmax=max_gpp_all, ax=axs[3, 0])
axs[3, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
plt.title('')
axs[3, 0].axis('off')
axs[3, 0].text(-0.07, 0.55, 'Clim$_{OSHD 0.5\degree}$+LU$_{Gl 0.5\degree}$', va='bottom', ha='center',
               rotation='vertical', rotation_mode='anchor', transform=axs[3, 0].transAxes, fontsize=17)
p1 = (gpp_oshd_05_avg.where((gpp_oshd_avg.mask_lake == 1) & (gpp_crujra_025_avg.mask_nan == 1) &
                            (gpp_crujra_025_avg.mask_nan_1km == 1)) -
     gpp_oshd_global_avg.where((gpp_oshd_avg.mask_lake == 1) &
                               (gpp_crujra_025_avg.mask_nan == 1) &
                               (gpp_crujra_025_avg.mask_nan_1km == 1))).plot(transform=proj_data, cmap='RdBu_r',
                                                                             add_colorbar=False, vmin=-diff_lim,
                                                                             vmax=diff_lim, ax=axs[3, 1])
axs[3, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
plt.title('')
plt.axis('off')
axs[3, 1].axis('off')

cbar_abs = inset_axes(axs[2, 0], width="70%",  # width: 5% of parent_bbox width
                      height="45%",  # height: 50%
                      loc="lower center", bbox_to_anchor=(0.5, -1.12, 1, 0.15), bbox_transform=axs[2, 0].transAxes,
                      borderpad=0.2)
cbar_delta = inset_axes(axs[2, 0], width="70%",  # width: 5% of parent_bbox width
                        height="45%",  # height: 50%
                        loc="lower center", bbox_to_anchor=(0.5, -1.45, 1, 0.15), bbox_transform=axs[2, 0].transAxes,
                        borderpad=0.2)
fig.colorbar(p, cax=cbar_abs, orientation="horizontal", label='GPP$_{JA}$ [gC m$^{-2}$ month$^{-1}$]')
fig.colorbar(p1, cax=cbar_delta, orientation="horizontal", label='$\Delta$ GPP$_{JA}$ [gC m$^{-2}$ month$^{-1}$]')
plt.subplots_adjust(right=0.99, top=0.99, wspace=-0.02, hspace=-0.33)
plt.tight_layout()
plt.show()
fig.savefig(Path(r'/home/lud11/malle/CLM5_CH/new_figures/gpp_spatial_comp_v2.pdf'),
            facecolor='white', transparent=False)
fig.savefig(Path(r'/home/lud11/malle/CLM5_CH/new_figures/gpp_spatial_comp_v2.png'),
            facecolor='white', transparent=False)
