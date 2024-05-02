# -*- coding: utf-8 -*-
"""
Desc:
Created on 08.12.22 18:02
@author: malle
"""

# run with python39_pickle environment!!!

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import cartopy.crs as ccrs
import cmocean
import matplotlib
import platform
import cartopy.io.shapereader as shapereader
from pyproj import Transformer
import pandas as pd
import seaborn as sns
import cftime

# set font type for PDFs and EPSs
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# make switch to windows if working from home
if platform.system() == 'Windows':
    base_dir = Path('L:\malle\CLM5_CH')
else:
    base_dir = Path('/home/lud11/malle/CLM5_CH')

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

highres_pai = input_highres.MONTHLY_LAI.mean(dim='lsmpft').sel(time=slice(0, 3)).mean(dim='time') + \
              input_highres.MONTHLY_SAI.mean(dim='lsmpft').sel(time=slice(0, 3)).mean(dim='time')
global_pai = input_global.MONTHLY_LAI.mean(dim='lsmpft').sel(time=slice(0, 3)).mean(dim='time') + \
             input_global.MONTHLY_SAI.mean(dim='lsmpft').sel(time=slice(0, 3)).mean(dim='time')

(highres_pai - global_pai).plot()
plt.show()

max_all_lai = np.max([input_highres.MONTHLY_LAI.max().values, input_global.MONTHLY_LAI.max().values])

# prep for plots
# Find the swiss boundary polygon.
countries = shapereader.natural_earth(resolution='10m', category='cultural', name='admin_0_countries')

for country in shapereader.Reader(countries).records():
    if 'Switzerland' in country.attributes['NAME_EN']:
        switzerland = country.geometry
        break
else:
    raise ValueError('Unable to find the CH boundary.')

proj_data = ccrs.PlateCarree()
proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)

fig1 = plt.figure(figsize=[8, 9])
ax0a = plt.subplot(1, 1, 1, projection=proj_map)
plot_pft = (highres_pai - global_pai)
p = plot_pft.plot(transform=proj_data, cmap='PiYG', center=0, add_colorbar=True, subplot_kws={'projection': proj_map},
                  cbar_kwargs=dict(location="right", label=r"$\Delta$PAI (LAI+SAI)", shrink=0.5))
ax0a.add_geometries([switzerland], ccrs.Geodetic(), edgecolor='darkred', linewidth=1.35, facecolor='none')
plt.axis('off')
ax0a.set_title('PAI difference between LU$_{HR1km}$ and LU$_{Gl1km}$', fontsize=11)
plt.tight_layout()
plt.show()
fig1.savefig(base_dir / 'figures_revisions' / 'delta_pai.png', facecolor='white', transparent=False, bbox_inches='tight')


# now load snow data
bf_fsm = base_dir / 'FSM_new/analysed_grid_python_out'
#FSM_2018 = xr.open_dataset(bf_fsm / 'analysis_2018_1km.nc')
FSM_int = xr.open_dataset(bf_fsm / 'snow_2018_hyps.nc')
FSM_feb = FSM_int.DATA.isel(time=1)

fig1 = plt.figure(figsize=[8, 9])
ax0a = plt.subplot(1, 1, 1, projection=proj_map)
p = FSM_feb.plot(transform=proj_data, cmap='viridis', add_colorbar=True, subplot_kws={'projection': proj_map},
                  cbar_kwargs=dict(location="right", label=r"$\Delta$PAI (LAI+SAI)", shrink=0.5))
ax0a.add_geometries([switzerland], ccrs.Geodetic(), edgecolor='darkred', linewidth=1.35, facecolor='none')
plt.axis('off')
ax0a.set_title('Snow depth as simulated from FSM2 for 1st Feb 2018', fontsize=11)
plt.tight_layout()
plt.show()
fig1.savefig(base_dir / 'figures_revisions' / 'snow_fsm2_feb2018.png', facecolor='white', transparent=False, bbox_inches='tight')

run_in = 'OSHD_FILES'
file_out = base_dir / run_in / 'HS_hyps.nc'
snow_in_sel = xr.open_dataset(file_out)
clm5_feb = snow_in_sel.DATA.sel(time=cftime.DatetimeNoLeap(2018, 2, 1, 0, 0, 0, 0, has_year_zero=True))

run_in = 'OSHD_FILES_OLD'
file_out = base_dir / run_in / 'HS_hyps.nc'
snow_in_sel = xr.open_dataset(file_out)
clm5_feb_gl = snow_in_sel.DATA.sel(time=cftime.DatetimeNoLeap(2018, 2, 1, 0, 0, 0, 0, has_year_zero=True))

fig1 = plt.figure(figsize=[8, 9])
ax0a = plt.subplot(1, 1, 1, projection=proj_map)
p = clm5_feb.plot(transform=proj_data, cmap='Purples', add_colorbar=True, subplot_kws={'projection': proj_map},
                  cbar_kwargs=dict(location="right", label=r"HS [m]", shrink=0.5))
ax0a.add_geometries([switzerland], ccrs.Geodetic(), edgecolor='darkred', linewidth=1.35, facecolor='none')
plt.axis('off')
ax0a.set_title('Snow depth as simulated by CLM5 Clim$_{OSHD1km}$LU$_{HR1km}$ for 1st Feb 2018', fontsize=11)
plt.tight_layout()
plt.show()
fig1.savefig(base_dir / 'figures_revisions' / 'snow_fsm2_feb2018.png', facecolor='white', transparent=False, bbox_inches='tight')


fig1 = plt.figure(figsize=[13, 7])
ax0a = plt.subplot(1, 2, 1, projection=proj_map)
p = plot_pft.plot(transform=proj_data, cmap='PiYG', center=0, add_colorbar=True, subplot_kws={'projection': proj_map},
                  cbar_kwargs=dict(location="right", label=r"$\Delta$PAI (LAI+SAI)", shrink=0.5))
ax0a.add_geometries([switzerland], ccrs.Geodetic(), edgecolor='darkred', linewidth=1.35, facecolor='none')
plt.axis('off')
ax0a.set_title('PAI difference between LU$_{HR1km}$ and LU$_{Gl1km}$', fontsize=11)

ax0a = plt.subplot(1, 2, 2, projection=proj_map)
p = clm5_feb.plot(transform=proj_data, cmap='Purples', add_colorbar=True, subplot_kws={'projection': proj_map},
                  cbar_kwargs=dict(location="right", label=r"HS [m]", shrink=0.5))
ax0a.add_geometries([switzerland], ccrs.Geodetic(), edgecolor='darkred', linewidth=1.35, facecolor='none')
plt.axis('off')
ax0a.set_title('Snow depth as simulated by CLM5 Clim$_{OSHD1km}$LU$_{HR1km}$ for 1st Feb 2018', fontsize=11)
plt.tight_layout()
plt.show()
fig1.savefig(base_dir / 'figures_revisions' / 'pai_snow_fsm2_feb2018.png', facecolor='white', transparent=False,
             bbox_inches='tight')

mask = (np.abs(plot_pft) > 0.25)

hr_large = clm5_feb.where(mask.data).data
hr_small = clm5_feb.where(~mask.data).data
gl_large = clm5_feb_gl.where(mask.data).data
gl_small = clm5_feb_gl.where(~mask.data).data

data_all = {'Clim$_{OSHD1km}$LU$_{HR1km}$, large $\Delta$LU': hr_large.flatten(),
        'Clim$_{OSHD1km}$LU$_{Gl1km}$, large $\Delta$LU': gl_large.flatten(),
     'Clim$_{OSHD1km}$LU$_{HR1km}$, small $\Delta$LU': hr_small.flatten(),
     'Clim$_{OSHD1km}$LU$_{Gl1km}$, small $\Delta$LU': gl_small.flatten()}

data_large = {'Clim$_{OSHD1km}$LU$_{HR1km}$, large $\Delta$LU': hr_large.flatten(),
        'Clim$_{OSHD1km}$LU$_{Gl1km}$, large $\Delta$LU': gl_large.flatten()}

data_small = {'Clim$_{OSHD1km}$LU$_{HR1km}$, small $\Delta$LU': hr_small.flatten(),
     'Clim$_{OSHD1km}$LU$_{Gl1km}$, small $\Delta$LU': gl_small.flatten()}

df_box = pd.DataFrame(data=data_all)
df_box_sm = pd.DataFrame(data=data_large)
df_box_lrg = pd.DataFrame(data=data_small)

fig = plt.figure()  # figsize=(12, 12)
flierprops = dict(marker='o', markeredgecolor='gray', markerfacecolor='silver', alpha=0.45)
axes = fig.add_subplot(111)
ax = sns.boxplot(data=df_box,
                 ax=axes,flierprops=flierprops, linewidth=1.6, saturation=0.8, showfliers=False,
                 meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black"},
                 medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), showmeans=True)
ax.set_ylabel('HS [m]')
plt.xticks(rotation=40, ha='right')
plt.tight_layout()
plt.show()
fig.savefig(base_dir / 'figures_revisions' / 'boxplot_comp.png', facecolor='white', transparent=False,
             bbox_inches='tight')

fig = plt.figure()  # figsize=(12, 12)
flierprops = dict(marker='o', markeredgecolor='gray', markerfacecolor='silver', alpha=0.45)
axes = fig.add_subplot(121)
ax = sns.boxplot(data=df_box_sm,
                 ax=axes,flierprops=flierprops, linewidth=1.6, saturation=0.8, showfliers=False,
                 meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black"},
                 medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), showmeans=True)
ax.set_ylabel('HS [m]')
plt.xticks(rotation=40, ha='right')

axes = fig.add_subplot(122)
ax = sns.boxplot(data=df_box_lrg,
                 ax=axes,flierprops=flierprops, linewidth=1.6, saturation=0.8, showfliers=False,
                 meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black"},
                 medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), showmeans=True)
ax.set_ylabel('HS [m]')
plt.xticks(rotation=40, ha='right')

plt.tight_layout()
plt.show()
fig.savefig(base_dir / 'figures_revisions' / 'boxplot_comp_sm_lrg.png', facecolor='white', transparent=False,
             bbox_inches='tight')

fig1 = plt.figure(figsize=[15, 8])
fig = plt.figure(figsize=(10, 7.5))
gs = fig.add_gridspec(3, 4)
ax1 = fig.add_subplot(gs[:2, :2], projection=proj_map)
ax2 = fig.add_subplot(gs[:2, 2:], projection=proj_map)
ax3 = fig.add_subplot(gs[2:, 1:3])

p = plot_pft.plot(ax=ax1, transform=proj_data, cmap='PiYG', center=0, add_colorbar=True,
                  cbar_kwargs=dict(location="right", label=r"$\Delta$PAI (LAI+SAI)", shrink=0.5))
ax1.add_geometries([switzerland], ccrs.Geodetic(), edgecolor='darkred', linewidth=1.35, facecolor='none')
ax1.axis('off')
ax1.set_title('(a)', loc='left')

#ax1.set_title('PAI difference between LU$_{HR1km}$ and LU$_{Gl1km}$', fontsize=11)

p = clm5_feb.plot(ax=ax2, transform=proj_data, cmap='Purples', add_colorbar=True,
                  cbar_kwargs=dict(location="right", label=r"HS [m]", shrink=0.5))
ax2.add_geometries([switzerland], ccrs.Geodetic(), edgecolor='darkred', linewidth=1.35, facecolor='none')
ax2.axis('off')
ax2.set_title('')
ax2.set_title('(b)', loc='left')
              #'Snow depth as simulated by CLM5 Clim$_{OSHD1km}$LU$_{HR1km}$ for 1st Feb 2018', fontsize=11)

ax = sns.boxplot(data=df_box,
                 ax=ax3,flierprops=flierprops, linewidth=1.6, saturation=0.8, showfliers=False,
                 meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black"},
                 medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), showmeans=True)
ax3.set_ylabel('HS [m]')
ax3.set_title('(c)', loc='left')
plt.xticks(rotation=20, ha='right')
plt.tight_layout()
plt.show()
fig.savefig(base_dir / 'figures_revisions' / 'comb_fig_lc.png', facecolor='white', transparent=False,
             bbox_inches='tight')