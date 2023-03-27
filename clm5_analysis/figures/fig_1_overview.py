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
import cmocean
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
veg_global_1km = input_global.PCT_NATVEG.where(domain.mask.data, np.nan)
veg_highres_1km = input_highres.PCT_NATVEG.where(domain.mask.data, np.nan)

bf_files = base_dir / 'out_for_plots'

# output data, so I can work with it on my windows machine
file_oshd_05deg = bf_files / 'temp_oshd_05deg.nc'
file_cru_05deg = bf_files / 'temp_cru_05deg.nc'
file_cruplus_05deg = bf_files / 'temp_cruplus_05deg.nc'
file_oshd_025deg = bf_files / 'temp_oshd_025deg.nc'
file_cru_025deg = bf_files / 'temp_cru_025deg.nc'
file_cruplus_025deg = bf_files / 'temp_cruplus_025deg.nc'
file_oshd_1km = bf_files / 'temp_oshd_1km.nc'
file_cru_1km = bf_files / 'temp_cru_1km.nc'
file_cruplus_1km = bf_files / 'temp_cruplus_1km.nc'

temp_oshd_05deg = xr.open_dataset(file_oshd_05deg)
temp_oshd_05deg['lon']=temp_oshd_05deg['lon']-180
temp_cru_05deg = xr.open_dataset(file_cru_05deg)
temp_cruplus_05deg = xr.open_dataset(file_cruplus_05deg)

temp_oshd_025deg = xr.open_dataset(file_oshd_025deg)
temp_oshd_025deg['lon'] = temp_oshd_025deg['lon']-180
temp_cru_025deg = xr.open_dataset(file_cru_025deg)
temp_cruplus_025deg = xr.open_dataset(file_cruplus_025deg)

temp_oshd_1km = xr.open_dataset(file_oshd_1km)
temp_cru_1km = xr.open_dataset(file_cru_1km)
temp_cruplus_1km = xr.open_dataset(file_oshd_1km)

# now vegetation
file_vegHR_025deg = bf_files / 'veg_gl_025.nc'
file_vegHR_05deg = bf_files / 'veg_gl_05.nc'
file_vegG_025deg = bf_files / 'veg_gl_025.nc'
file_vegG_05deg = bf_files / 'veg_gl_05.nc'

veg_highres_025deg = xr.open_dataset(file_vegHR_025deg)
veg_highres_05deg = xr.open_dataset(file_vegHR_05deg)
veg_global_025deg = xr.open_dataset(file_vegG_025deg)
veg_global_05deg = xr.open_dataset(file_vegG_05deg)

test_oshd = xr.open_dataset(bf_files / 'test_oshd_05deg.nc')

base_dir_t = Path('C:\\Users\malle\Documents')
input_oshd = xr.open_dataset(base_dir_t / 'clmforc.OSHD1km.TQ.2018-05.nc')
month_avg_oshd = input_oshd.mean(dim='time')
del input_oshd
input_crujra = xr.open_dataset(base_dir_t / 'clmforc.crujra.TQ.2018-05.nc')
month_avg_cru = input_crujra.mean(dim='time')
del input_crujra

# now plotting etc.
# plotting settings
proj_data = ccrs.PlateCarree()
proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)
fig1 = plt.figure(figsize=[9, 6.5])
ax0a = plt.subplot(111, projection=proj_map)
#veg_highres_025deg['__xarray_dataarray_variable__'].plot(transform=proj_data, cmap=cmocean.cm.thermal, subplot_kws={'projection': proj_map})
temp_oshd_025deg.TBOT.plot(transform=proj_data, cmap='viridis', subplot_kws={'projection': proj_map})
ax0a.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
plt.axis('off')
plt.show()

# start with comp of CRUJRA data
proj_data = ccrs.PlateCarree()
proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)

min_all = np.nanmin([temp_cru.TBOT-273.15, temp_month_avg_cru_025deg-273.15, temp_month_avg_cru_1deg-273.15, temp_cru_plus.TBOT-273.15, temp_month_avg_cruplus_025deg-273.15, temp_month_avg_cruplus_1deg-273.15,temp_oshd.TBOT-273.15])
max_all = np.nanmax([temp_cru.TBOT-273.15, temp_month_avg_cru_025deg-273.15, temp_month_avg_cru_1deg-273.15, temp_cru_plus.TBOT-273.15, temp_month_avg_cruplus_025deg-273.15, temp_month_avg_cruplus_1deg-273.15,temp_oshd.TBOT-273.15])

fig1 = plt.figure(figsize=[9, 6.5])
ax0a = plt.subplot(331, projection=proj_map)
tbot_c = temp_month_avg_cru_1deg - 273.15
p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
                vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
ax0a.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
plt.axis('off')

ax1 = plt.subplot(332, projection=proj_map)
tbot_c = temp_month_avg_cru_025deg - 273.15
p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
                vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
ax1.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
plt.axis('off')

ax0 = plt.subplot(333, projection=proj_map)
tbot_c = temp_cru.TBOT - 273.15
p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
                vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
ax0.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
plt.axis('off')

ax0a = plt.subplot(334, projection=proj_map)
tbot_c = temp_month_avg_cruplus_1deg - 273.15
p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
                vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
ax0a.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
plt.axis('off')

ax1 = plt.subplot(335, projection=proj_map)
tbot_c = temp_month_avg_cruplus_025deg - 273.15
p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
                vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
ax1.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
plt.axis('off')

ax0 = plt.subplot(336, projection=proj_map)
tbot_c = temp_cru_plus.TBOT - 273.15
p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
                vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
ax0.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
plt.axis('off')

cax = fig1.add_axes([ax0.get_position().x1+0.01, ax0.get_position().y0-0.08, 0.022, ax0.get_position().height*2])
cb = plt.colorbar(p, cax=cax)
cb.ax.tick_params(labelsize=8)
cb.set_label('Air Temperature ($\degree$C)', rotation=90, fontsize=9)

ax0a = plt.subplot(337, projection=proj_map)
tbot_c = temp_month_avg_oshd_1deg - 273.15
p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
                vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
ax0a.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=0.9)
plt.axis('off')

ax0a = plt.subplot(338, projection=proj_map)
tbot_c = temp_month_avg_oshd_025deg - 273.15
p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
                vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
ax0a.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=0.9)
plt.axis('off')

ax0a = plt.subplot(339, projection=proj_map)
tbot_c = temp_oshd.TBOT - 273.15
p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
                vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
ax0a.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=0.9)
plt.axis('off')
plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()

# fig1.savefig(base_dir / "test.svg", bbox_inches='tight', format='svg')
fig1.savefig(base_dir / 'comp_temp_cru_all.pdf')
fig1.savefig(base_dir / 'comp_temp_cru_all.png')

fig1.savefig(base_dir / 'comp_temp_cru_all.eps')

plt.show()
#
#
# # now surface dataset
# min_all = 0
# max_all = 100
#
# fig1 = plt.figure(figsize=[9, 6.5])
# ax0a = plt.subplot(331, projection=proj_map)
# tbot_c = veg_global_1deg_back1km_dom
# p = tbot_c.plot(transform=proj_data, cmap='viridis', add_colorbar=False,
#                 vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
# ax0a.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
# plt.axis('off')
#
# ax1 = plt.subplot(332, projection=proj_map)
# tbot_c = veg_global_025deg_back1km_dom
# p = tbot_c.plot(transform=proj_data, cmap='viridis', add_colorbar=False,
#                 vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
# ax1.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
# plt.axis('off')
#
# ax0 = plt.subplot(333, projection=proj_map)
# tbot_c = veg_global
# p = tbot_c.plot(transform=proj_data, cmap='viridis', add_colorbar=False,
#                 vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
# ax0.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
# plt.axis('off')
#
# ax0a = plt.subplot(334, projection=proj_map)
# tbot_c = veg_highres_05deg
# p = tbot_c.plot(transform=proj_data, cmap='viridis', add_colorbar=False,
#                 vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
# ax0a.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
# plt.axis('off')
#
# ax0a = plt.subplot(335, projection=proj_map)
# tbot_c = veg_highres_025deg
# p = tbot_c.plot(transform=proj_data, cmap='viridis', add_colorbar=False,
#                 vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
# ax0a.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
# plt.axis('off')
#
# ax0a = plt.subplot(336, projection=proj_map)
# tbot_c = veg_highres
# p = tbot_c.plot(transform=proj_data, cmap='viridis', add_colorbar=False,
#                 vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
# ax0a.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
# plt.axis('off')
#
# cax = fig1.add_axes([ax0.get_position().x1+0.01, ax0.get_position().y0-0.22, 0.02, ax0.get_position().height*1.5])
# cb = plt.colorbar(p, cax=cax)
# cb.ax.tick_params(labelsize=8)
# cb.set_label('Vegetation Cover (%)', rotation=90, fontsize=9)
#
# plt.tight_layout()
# plt.subplots_adjust(wspace=0, hspace=0)
# plt.show()
# # fig1.savefig(base_dir / "test.svg", bbox_inches='tight', format='svg')
# fig1.savefig(base_dir / 'comp_vcover_cru_all.pdf')
# fig1.savefig(base_dir / 'comp_vcover_cru_all.png')
#
# fig1.savefig(base_dir / 'comp_vcover_cru_all.eps')
#
# # use this if only CH:
# #proj_map = ccrs.epsg(2056)
# proj_data = ccrs.PlateCarree()
# proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)
#
# min_all = np.nanmin([temp_cru.TBOT-273.15, temp_oshd.TBOT-273.15, temp_cru_plus.TBOT-273.15])
# max_all = np.nanmax([temp_cru.TBOT-273.15, temp_oshd.TBOT-273.15, temp_cru_plus.TBOT-273.15])
#
# fig1 = plt.figure(figsize=[8, 1.5])
# ax0a = plt.subplot(141, projection=proj_map)
# tbot_c = temp_cru.TBOT - 273.15
# p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
#                 vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
# ax0a.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
# plt.axis('off')
#
# ax1 = plt.subplot(142, projection=proj_map)
# tbot_c = temp_cru_plus.TBOT - 273.15
# p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
#                 vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
# ax1.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
# plt.axis('off')
#
# ax0 = plt.subplot(143, projection=proj_map)
# tbot_c = temp_oshd.TBOT - 273.15
# p = tbot_c.plot(transform=proj_data, cmap=cmocean.cm.thermal, add_colorbar=False,
#                 vmin=min_all, vmax=max_all, subplot_kws={'projection': proj_map})
# ax0.add_feature(cf.BORDERS, linewidth=1.2, edgecolor='dimgray', alpha=1)
# plt.axis('off')
#
# plt.tight_layout()
# cax = fig1.add_axes([ax0.get_position().x1+0.04, ax0.get_position().y0, 0.02, ax0.get_position().height])
# cb = plt.colorbar(p, cax=cax)
# cb.ax.tick_params(labelsize=8)
# cb.set_label('Air Temperature ($\degree$C)', rotation=90, fontsize=9)
# plt.show()
# # fig1.savefig(base_dir / "test.svg", bbox_inches='tight', format='svg')
# fig1.savefig(Path(r'C:\Users\malle\Documents\paper_clm5_ch\temp_comp.pdf'), transparent=True)
# plt.show()
#
#
# transformer = Transformer.from_crs("epsg:4326", "epsg:32632")
# x2, y2 = transformer.transform(month_avg_cruplus.LATIXY,month_avg_cruplus.LONGXY)
#
# ccrs.UTM(zone=32, southern_hemisphere=False)
# proj_map = ccrs.epsg(32632)
#
# # now make grid
# stations_1000 = pd.read_csv(base_dir / 'coor_xy_1000.xy', delimiter=' ', header=None, names=['x', 'y'])
# stations_2000 = pd.read_csv(base_dir / 'coor_xy_2000.xy', delimiter=' ', header=None, names=['x', 'y'])
# stations_3000 = pd.read_csv(base_dir / 'coor_xy_3000.xy', delimiter=' ', header=None, names=['x', 'y'])
# fluxnet_pts = pd.read_csv(base_dir / 'fluxnet.xy', delimiter=' ', header=None, names=['x', 'y'])
# fig = plt.figure(figsize=[9, 6.5])
# ax = plt.subplot(331, projection=proj_map)
# major_ticks_x = np.linspace(np.min(x2), np.max(x2),11)
# major_ticks_y = np.linspace(np.min(y2), np.max(y2),7)
# ax.set_xticks(major_ticks_x)
# ax.set_yticks(major_ticks_y)
# ax.set_axisbelow(True)
# ax.grid(which='both')
# ax.grid(which='major', linewidth=0.5, color='darkred', alpha=0.7)
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# ax.tick_params(direction='in', length=0, width=0)
# ax.add_feature(cf.BORDERS, linewidth=0.9, edgecolor='black', alpha=1)
# ax.set_extent([np.min(lon_geo), np.max(lon_geo), np.min(lat_geo),np.max(lat_geo)])
# #plt.scatter(stations_1000.x, stations_1000.y, marker='.', s=13, color='skyblue', transform=proj_data, edgecolor='black', linewidth=0.18)
# #plt.scatter(stations_2000.x, stations_2000.y, marker='.', s=13, color='skyblue', transform=proj_data, edgecolor='black', linewidth=0.18)
# #plt.scatter(stations_3000.x, stations_3000.y, marker='.', s=13, color='skyblue', transform=proj_data, edgecolor='black', linewidth=0.18)
# plt.scatter(fluxnet_pts.x, fluxnet_pts.y, marker='.', s=14, color='forestgreen', transform=proj_data, edgecolor='black', linewidth=0.15)
#
# ax = fig.add_subplot(332, projection=proj_map)
# major_ticks_x = np.linspace(np.min(x2), np.max(x2), 20)
# major_ticks_y = np.linspace(np.min(y2), np.max(y2), 12)
# ax.set_xticks(major_ticks_x)
# ax.set_yticks(major_ticks_y)
# ax.set_axisbelow(True)
# ax.grid(which='both')
# ax.grid(which='major', linewidth=0.4, color='darkred', alpha=0.7)
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# ax.tick_params(direction='in', length=0, width=0)
# ax.add_feature(cf.BORDERS, linewidth=0.9, edgecolor='black', alpha=1)
# ax.set_extent([np.min(lon_geo), np.max(lon_geo), np.min(lat_geo), np.max(lat_geo)])
# #plt.scatter(fluxnet_pts.x, fluxnet_pts.y, marker='.', s=13, color='forestgreen', transform=proj_data, edgecolor='black', linewidth=0.18)
#
# ax = fig.add_subplot(333, projection=proj_map)
# major_ticks_x = np.linspace(np.min(x2), np.max(x2),366)
# major_ticks_y = np.linspace(np.min(y2),np.max(y2),273)
# ax.set_xticks(major_ticks_x)
# ax.set_yticks(major_ticks_y)
# ax.set_axisbelow(True)
# ax.grid(which='both')
# ax.grid(which='major', linewidth=0.2, color='darkred', alpha=0.7)
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# ax.tick_params(direction='in', length=0, width=0)
# ax.add_feature(cf.BORDERS, linewidth=0.9, edgecolor='black', alpha=1)
# ax.set_extent([np.min(lon_geo), np.max(lon_geo), np.min(lat_geo),np.max(lat_geo)])
#
# plt.scatter(stations_1000.x, stations_1000.y, marker='.', s=13, color='skyblue', transform=proj_data, edgecolor='black', linewidth=0.15)
# plt.scatter(stations_2000.x, stations_2000.y, marker='.', s=13, color='skyblue', transform=proj_data, edgecolor='black', linewidth=0.15)
# plt.scatter(stations_3000.x, stations_3000.y, marker='.', s=13, color='skyblue', transform=proj_data, edgecolor='black', linewidth=0.15)
# plt.scatter(fluxnet_pts.x, fluxnet_pts.y, marker='.', s=13, color='forestgreen', transform=proj_data, edgecolor='black', linewidth=0.15)
#
# ax.add_feature(cf.BORDERS, linewidth=0.9, edgecolor='black', alpha=1)
# ax.set_extent([np.min(lon_geo), np.max(lon_geo), np.min(lat_geo), np.max(lat_geo)])
#
# plt.show()
# fig.savefig(base_dir / 'grid_comp2.pdf', transparent=True)
# fig.savefig(base_dir / 'grid_comp2.png')
# fig.savefig(base_dir / 'grid_comp2.eps')
#
#
# fig = plt.figure(figsize=[8, 1.5])
# ax = fig.add_subplot(141, projection=proj_map)
# major_ticks_x = np.linspace(np.min(x2),np.max(x2),7)
# major_ticks_y = np.linspace(np.min(y2),np.max(y2),4)
# ax.set_xticks(major_ticks_x)
# ax.set_yticks(major_ticks_y)
# ax.set_axisbelow(True)
# ax.grid(which='both')
# ax.grid(which='major',linewidth=0.5, color='darkred', alpha=0.7)
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# ax.tick_params(direction='in', length=0, width=0)
# ax.add_feature(cf.BORDERS, linewidth=0.9, edgecolor='black', alpha=1)
# ax.set_extent([np.min(lon_geo),np.max(lon_geo),np.min(lat_geo),np.max(lat_geo)])
#
# ax = fig.add_subplot(142, projection=proj_map)
# major_ticks_x = np.linspace(np.min(x2),np.max(x2),20)
# major_ticks_y = np.linspace(np.min(y2),np.max(y2),11)
# ax.set_xticks(major_ticks_x)
# ax.set_yticks(major_ticks_y)
# ax.set_axisbelow(True)
# ax.grid(which='both')
# ax.grid(which='major',linewidth=0.4, color='darkred', alpha=0.7)
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# ax.tick_params(direction='in', length=0, width=0)
# ax.add_feature(cf.BORDERS, linewidth=0.9, edgecolor='black', alpha=1)
# ax.set_extent([np.min(lon_geo),np.max(lon_geo),np.min(lat_geo),np.max(lat_geo)])
#
# ax = fig.add_subplot(143, projection=proj_map)
# major_ticks_x = np.linspace(np.min(x2),np.max(x2),366)
# major_ticks_y = np.linspace(np.min(y2),np.max(y2),273)
# ax.set_xticks(major_ticks_x)
# ax.set_yticks(major_ticks_y)
# ax.set_axisbelow(True)
# ax.grid(which='both')
# ax.grid(which='major',linewidth=0.2, color='darkred', alpha=0.7)
# ax.set_xticklabels([])
# ax.set_yticklabels([])
# ax.tick_params(direction='in', length=0, width=0)
# ax.add_feature(cf.BORDERS, linewidth=0.9, edgecolor='black', alpha=1)
# ax.set_extent([np.min(lon_geo),np.max(lon_geo),np.min(lat_geo),np.max(lat_geo)])
#
# ax = fig.add_subplot(144, projection=proj_map)
# plt.scatter(stations_1000.x, stations_1000.y, marker='.', s=9, color='darkred', transform=proj_data)
# plt.scatter(stations_2000.x, stations_2000.y, marker='.', s=9, color='darkred', transform=proj_data)
# plt.scatter(stations_3000.x, stations_3000.y, marker='.', s=9, color='darkred', transform=proj_data)
# ax.add_feature(cf.BORDERS, linewidth=0.9, edgecolor='black', alpha=1)
# ax.set_extent([np.min(lon_geo), np.max(lon_geo), np.min(lat_geo), np.max(lat_geo)])
# #plt.axis('off')
#
# plt.show()
# fig.savefig(Path(r'C:\Users\malle\Documents\paper_clm5_ch\grid_comp.pdf'), transparent=True)
#
#
#
# fig = plt.figure(figsize=[8, 1.5])
#
#
# plt.show()
#
# gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=2, color='gray', alpha=0.5, linestyle='--')
# gl.xlabels_top = False
# gl.ylabels_left = False
# gl.xlocator = mticker.FixedLocator(major_ticks_x)
# gl.ylocator = mticker.FixedLocator(major_ticks_y)
# #gl.xlines = False
# #gl.ylines = False
# gl.xlocator = mticker.FixedLocator(major_ticks_x)
# gl.ylocator = mticker.FixedLocator(major_ticks_y)
#
# # if grid lines wanted this could be useful:
# # gl = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
# #                   linewidth=2, color='gray', alpha=0.5, linestyle='--')
# # gl.xlabels_top = False
# # gl.ylabels_left = False
# # gl.xlines = False
# # gl.ylines = False
# # gl.xlocator = mticker.FixedLocator([7, 8, 9, 10])
# # gl.ylocator = mticker.FixedLocator([46,46.5, 47, 47.5])
# #
# # gl.xformatter = LONGITUDE_FORMATTER
# # gl.yformatter = LATITUDE_FORMATTER
# # gl.xlabel_style = {'size': 15, 'color': 'gray'}
# # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
#
# #another way of transforming data:
#
# inProj = Proj('epsg:4326')
# outProj = Proj('epsg:32632')
# transformer = Transformer.from_crs("epsg:4326", "epsg:32632")
# x2, y2 = transformer.transform(month_avg_cruplus.LATIXY,month_avg_cruplus.LONGXY)
# #
# #
# # fig=plt.figure(figsize=(6, 3))
# # ax = plt.subplot(111, projection=proj_map)
# # ax.plot(y2, x2, temp_month_avg_cruplus.data)
# # #ax.add_feature(cf.BORDERS)
# # plt.show()