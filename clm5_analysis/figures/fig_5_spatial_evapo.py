# -*- coding: utf-8 -*-
"""
Desc:
Created on 13.01.23 16:01
@author: malle
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib
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

# create target grids for xesmf
domain = xr.open_dataset(bf / 'domain.lnd.CH_1km_navy.210407_new.nc')
surf_025deg = xr.open_dataset(bf / 'surfdata_0.25deg_CH_v2_hist_16pfts_Irrig_CMIP6_simyr1850_c210415.nc')
lat_new = surf_025deg.LATIXY.data[:, 0]
lon_new = surf_025deg.LONGXY.data[0, :]
ds_target_025 = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

surf_05deg = xr.open_dataset(bf / 'surfdata_CH_05deg_v2_hist_16pfts_Irrig_CMIP6_simyr2000_c230113.nc')
lat_new = surf_05deg.LATIXY.data[:, 0]
lon_new = surf_05deg.LONGXY.data[0, :]
ds_target_05 = xr.Dataset({"lat": (["lat"], lat_new), "lon": (["lon"], lon_new)})

# load clm5 simulations for the various model configurations -> focus on 2017 sum
evapo_oshd_avg = xr.open_dataset(bf / 'OSHD_FILES' / 'QFLX_EVAP_TOT_SP_YR_SUM.nc')
evapo_oshd_2017 = evapo_oshd_avg.isel(year=7)

evapo_oshd_global_avg = xr.open_dataset(bf / 'OSHD_FILES_OLD' / 'QFLX_EVAP_TOT_SP_YR_SUM.nc')
evapo_oshd_global_2017 = evapo_oshd_global_avg.isel(year=7)

evapo_crujra_avg = xr.open_dataset(bf / 'CRUJRA_FILES_OLD' / 'QFLX_EVAP_TOT_SP_YR_SUM.nc')
evapo_crujra_2017 = evapo_crujra_avg.isel(year=7)

evapo_crujra_nolapse_avg = xr.open_dataset(bf / 'CRUJRA_FILES_noLapse_OLD' / 'QFLX_EVAP_TOT_SP_YR_SUM.nc')
evapo_crujra_nolapse_2017 = evapo_crujra_nolapse_avg.isel(year=7)

evapo_crujra_025 = xr.open_dataset(bf / 'CRUJRA_FILES_025deg' / 'QFLX_EVAP_TOT_SP_YR_SUM.nc')
evapo_crujra_025_2017 = evapo_crujra_025.isel(year=7)

crujra_coarse = xr.open_dataset(bf / 'CRUJRA_FILES_05deg_cru_new' / 'QFLX_EVAP_TOT.nc') * 3600 * 24  # mm/s to mm/day
sum_yr = crujra_coarse.resample(time='1Y').sum()
bacK_regridder_nearest_05 = xe.Regridder(ds_target_05, evapo_crujra_nolapse_2017, "nearest_s2d")
sum_month_re = bacK_regridder_nearest_05(sum_yr)
evapo_crujra_05_2017 = sum_month_re.isel(time=7)

evapo_oshd_05 = xr.open_dataset(bf / 'OSHD_FILES_05_new' / 'QFLX_EVAP_TOT.nc')*3600*24  # mm/s to mm/day
sum_yr = evapo_oshd_05.resample(time='1Y').sum()
sum_month_re = bacK_regridder_nearest_05(sum_yr)
evapo_oshd_05_2017 = sum_month_re.isel(time=7)

bacK_regridder_nearest_025 = xe.Regridder(ds_target_025, evapo_oshd_2017, "nearest_s2d")
evapo_cr_025_2017_re = bacK_regridder_nearest_025(evapo_crujra_025_2017)

id_mask_1km = ~np.isnan(evapo_oshd_2017.evapo)
evapo_cr_025_2017_re.coords['mask_nan_1km'] = (('lat', 'lon'), (id_mask_1km*1).squeeze().data)

a = np.array([[evapo_oshd_2017.evapo], [evapo_oshd_global_2017.evapo], [evapo_crujra_2017.evapo],
              [evapo_crujra_nolapse_2017.evapo], [evapo_oshd_05_2017.DATA]])
a1 = np.array([[evapo_oshd_global_2017.evapo - evapo_oshd_2017.evapo],
               [evapo_crujra_2017.evapo - evapo_oshd_2017.evapo],
               [evapo_crujra_nolapse_2017.evapo - evapo_oshd_2017.evapo],
               [evapo_oshd_05_2017.DATA - evapo_oshd_2017.evapo]])

min_all, max_all = np.nanmin(a), np.nanmax(a)
min_diff, max_diff = np.nanmin(a1), np.nanmax(a1)
diff_lim = np.max([[np.abs(min_diff)], [np.abs(max_diff)]])

# now plot:
proj_data = ccrs.PlateCarree()
proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)

font = {'size': 15}
plt.rc('font', **font)

fig, axs = plt.subplots(4, 2, figsize=[9, 15], frameon=False, subplot_kw={'projection': proj_map})
evapo_oshd_global_2017.evapo.plot(transform=proj_data, cmap='YlGnBu', add_colorbar=False,
                                  vmin=min_all, vmax=max_all, ax=axs[0, 1])
axs[0, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
axs[0, 1].title.set_text('')
axs[0, 1].text(-0.07, 0.55, 'Clim$_{OSHD 1km}$+LU$_{Gl 1km}$', va='bottom', ha='center', rotation='vertical',
               rotation_mode='anchor', transform=axs[0, 1].transAxes, fontsize=17)
axs[0, 1].axis('off')
axs[0, 0].axis('off')
# effect of land use info:
evapo_oshd_2017.evapo.plot(transform=proj_data, cmap='YlGnBu', add_colorbar=False,
                           vmin=min_all, vmax=max_all, ax=axs[1, 0])
axs[1, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
axs[1, 0].title.set_text('')
axs[1, 0].axis('off')
axs[1, 0].text(-0.07, 0.55, 'Clim$_{OSHD 1km}$+LU$_{HR 1km}$', va='bottom', ha='center', rotation='vertical',
               rotation_mode='anchor', transform=axs[1, 0].transAxes, fontsize=17)
(evapo_oshd_2017.evapo - evapo_oshd_global_2017.evapo).\
     plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False, vmin=-diff_lim, vmax=diff_lim, ax=axs[1, 1])
axs[1, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
axs[1, 1].title.set_text('')
axs[1, 1].axis('off')
# effect of climate forcing:
evapo_crujra_nolapse_2017.evapo.where(evapo_cr_025_2017_re.mask_nan_1km == 1).plot(
     transform=proj_data, cmap='YlGnBu', add_colorbar=False, vmin=min_all, vmax=max_all, ax=axs[2, 0])
axs[2, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
axs[2, 0].title.set_text('')
axs[2, 0].axis('off')
axs[2, 0].text(-0.07, 0.55, 'Clim$_{CRU 1km}$+LU$_{Gl 1km}$', va='bottom', ha='center', rotation='vertical',
               rotation_mode='anchor', transform=axs[2, 0].transAxes, fontsize=17)

(evapo_crujra_nolapse_2017.evapo.where(evapo_cr_025_2017_re.mask_nan_1km == 1) -
 evapo_oshd_global_2017.evapo.where(evapo_cr_025_2017_re.mask_nan_1km == 1)).plot(transform=proj_data,
                                                                                  cmap='RdBu_r',
                                                                                  add_colorbar=False,
                                                                                  vmin=-diff_lim,
                                                                                  vmax=diff_lim,
                                                                                  ax=axs[2, 1])
axs[2, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
axs[2, 1].title.set_text('')
plt.axis('off')
axs[2, 1].axis('off')
# effect of resolution:
p = evapo_oshd_05_2017.DATA.where(evapo_cr_025_2017_re.mask_nan_1km == 1).plot(
     transform=proj_data, cmap='YlGnBu', add_colorbar=False, vmin=min_all, vmax=max_all, ax=axs[3, 0])
axs[3, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
axs[3, 0].title.set_text('')
axs[3, 0].axis('off')
axs[3, 0].text(-0.07, 0.55, 'Clim$_{OSHD 0.5\degree}$+LU$_{Gl 0.5\degree}$', va='bottom', ha='center',
               rotation='vertical', rotation_mode='anchor', transform=axs[3, 0].transAxes, fontsize=17)
p1 = (evapo_oshd_05_2017.DATA.where(evapo_cr_025_2017_re.mask_nan_1km == 1) -
      evapo_oshd_global_2017.evapo.where(evapo_cr_025_2017_re.mask_nan_1km == 1)).plot(transform=proj_data,
                                                                                       cmap='RdBu_r',
                                                                                       add_colorbar=False,
                                                                                       vmin=-diff_lim,
                                                                                       vmax=diff_lim,
                                                                                       ax=axs[3, 1])
axs[3, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
axs[3, 1].title.set_text('')
plt.axis('off')
axs[3, 1].axis('off')

cbar_abs = inset_axes(axs[2, 0], width="70%",  # width: 5% of parent_bbox width
                      height="45%",  # height: 50%
                      loc="lower center", bbox_to_anchor=(0.5, -1.12, 1, 0.15),
                      bbox_transform=axs[2, 0].transAxes, borderpad=0.2)
cbar_delta = inset_axes(axs[2, 0], width="70%",  # width: 5% of parent_bbox width
                        height="45%",  # height: 50%
                        loc="lower center", bbox_to_anchor=(0.5, -1.45, 1, 0.15),
                        bbox_transform=axs[2, 0].transAxes, borderpad=0.2)
fig.colorbar(p, cax=cbar_abs, orientation="horizontal", label='Evapotranspiration [mm]')
fig.colorbar(p1, cax=cbar_delta, orientation="horizontal", label='$\Delta$ Evapotranspiration [mm]')
plt.subplots_adjust(right=0.99, top=0.99, wspace=-0.02, hspace=-0.33)
plt.tight_layout()

fig.savefig(Path(r'/home/lud11/malle/CLM5_CH/new_figures/evapo_spatial_comp_v2.pdf'),
            facecolor='white', transparent=False)
fig.savefig(Path(r'/home/lud11/malle/CLM5_CH/new_figures/evapo_spatial_comp_v2.png'),
            facecolor='white', transparent=False)
plt.show()

# for supp. material: look at total precipiation
bf_precip = Path('/home/malle/Documents/input_checks/precip')
precip_2017_oshd = xr.open_dataset(bf_precip / 'yr_sum_OSHD_2017.nc')
precip_2017_cru = xr.open_dataset(bf_precip / 'yr_sum_CRUJRA_2017.nc')

precip_2017_oshd1 = precip_2017_oshd.isel(year=0).where(~np.isnan(evapo_oshd_2017.evapo.data), np.nan)
precip_2017_cru1 = precip_2017_cru.isel(year=0).where(~np.isnan(evapo_oshd_2017.evapo.data), np.nan)

precip_2017_cru1a = precip_2017_cru1.assign_coords({'lat': evapo_oshd_2017.lat})
precip_2017_cru1b = precip_2017_cru1a.assign_coords({'lon': evapo_oshd_2017.lon})
precip_2017_oshd1a = precip_2017_oshd1.assign_coords({'lat': evapo_oshd_2017.lat})
precip_2017_oshd1b = precip_2017_oshd1a.assign_coords({'lon': evapo_oshd_2017.lon})

# start with precip comp.
a = np.array([[precip_2017_oshd1.precip], [precip_2017_cru1.precip]])
min_all, max_all = np.nanmin(a), np.nanmax(a)

fig, axs = plt.subplots(1, 3, figsize=[8, 3], frameon=False, subplot_kw={'projection': proj_map})
p = precip_2017_oshd1b.precip.plot(transform=proj_data, cmap='YlGnBu', add_colorbar=False,
                                   vmin=min_all, vmax=max_all, ax=axs[0])
axs[0].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
axs[0].set_title('')
axs[0].axis('off')
axs[0].text(-0.05, 0.55, 'Clim$_{OSHD 1km}$', va='bottom', ha='center', rotation='vertical', rotation_mode='anchor',
            transform=axs[0].transAxes)
p1 = precip_2017_cru1b.precip.plot(transform=proj_data, cmap='YlGnBu',
                                   add_colorbar=False, vmin=min_all,
                                   vmax=max_all, ax=axs[1])
axs[1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
axs[1].set_title('')
axs[1].axis('off')
axs[1].text(-0.001, 0.55, 'Clim$_{CRU 1km}$', va='bottom', ha='center', rotation='vertical', rotation_mode='anchor',
            transform=axs[1].transAxes)

p2 = (precip_2017_cru1b.precip - precip_2017_oshd1b.precip).plot(transform=proj_data, cmap='RdBu_r',
                                                                 add_colorbar=False, ax=axs[2])
axs[2].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
axs[2].set_title('')
axs[2].axis('off')
axs[2].text(-0.001, 0.55, '$\Delta$ Clim$_{CRU 1km}$ - Clim$_{OSHD 1km}$', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',
            transform=axs[2].transAxes)

fig.colorbar(p1, ax=[axs[:2]], shrink=0.75, location='bottom', label='Yearly Precipitation 2017 [mm]')
fig.colorbar(p2, ax=[axs[2]], shrink=0.75, location='bottom', label='$\Delta$ Yearly Precipitation 2017 [mm]')

plt.tight_layout()
fig.savefig(Path(r'/home/lud11/malle/CLM5_CH/new_figures/precip_comp.png'))
fig.savefig(Path(r'/home/lud11/malle/CLM5_CH/new_figures/precip_comp.pdf'), transparent=True)
