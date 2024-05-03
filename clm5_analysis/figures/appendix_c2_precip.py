# -*- coding: utf-8 -*-
"""
Desc: 
Created on 03.08.23 09:40
@author: malle
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cf
import platform


# make switch to windows if working from laptop
if platform.system() == 'Windows':
    bf = Path('L:\malle\CLM5_CH')
else:
    bf = Path('/home/lud11/malle/CLM5_CH')

################################ Figure S8 #################################
# compare total precipiation for 2017 between forcing datasets OSHD & CRU

bf_precip = Path('/home/malle/Documents/input_checks/precip')

proj_data = ccrs.PlateCarree()
proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)

evapo_oshd_avg = xr.open_dataset(bf / 'OSHD_FILES' / 'QFLX_EVAP_TOT_SP_YR_SUM.nc')
evapo_oshd_2017 = evapo_oshd_avg.isel(year=7)
precip_2017_oshd = xr.open_dataset(bf_precip / 'yr_sum_OSHD_2017.nc')
precip_2017_cru = xr.open_dataset(bf_precip / 'yr_sum_CRUJRA_2017.nc')

precip_2017_oshd1 = precip_2017_oshd.isel(year=0).where(~np.isnan(evapo_oshd_2017.evapo.data), np.nan)
precip_2017_cru1 = precip_2017_cru.isel(year=0).where(~np.isnan(evapo_oshd_2017.evapo.data), np.nan)
precip_2017_cru1a = precip_2017_cru1.assign_coords({'lat': evapo_oshd_2017.lat})
precip_2017_cru1b = precip_2017_cru1a.assign_coords({'lon': evapo_oshd_2017.lon})
precip_2017_oshd1a = precip_2017_oshd1.assign_coords({'lat': evapo_oshd_2017.lat})
precip_2017_oshd1b = precip_2017_oshd1a.assign_coords({'lon': evapo_oshd_2017.lon})

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
fig.savefig(bf / 'new_figures' / 'precip_comp.png')
fig.savefig(bf / 'new_figures' / 'precip_comp.pdf', transparent=True)
