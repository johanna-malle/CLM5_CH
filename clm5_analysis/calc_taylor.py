# -*- coding: utf-8 -*-
"""
Desc: 
Created on 30.01.24 10:27
@author: malle
"""

from pathlib import Path
import xarray as xr
import glob
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import skill_metrics as sm
import numpy as np
import platform
import distinctipy
import cartopy.crs as ccrs
import cartopy.feature as cf

proj_data = ccrs.PlateCarree()
proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)


# make switch to windows if working from home
if platform.system() == 'Linux':
    bf = Path('/home/lud11/malle/CLM5_CH')
else:
    bf = Path('L:\malle\CLM5_CH')

bf_oshd = bf / 'FSM_new/analysed_grid_python_out'
all_FSM = xr.open_mfdataset(glob.glob(str(bf_oshd / '*025*')))

bf_out = bf / 'figures_revisions'


FSM_feb = all_FSM.isel(time=(all_FSM.time.dt.month == 2) & (all_FSM.time.dt.day == 1))
FSM_apr = all_FSM.isel(time=(all_FSM.time.dt.month == 4) & (all_FSM.time.dt.day == 1))
FSM_dez = all_FSM.isel(time=(all_FSM.time.dt.month == 12) & (all_FSM.time.dt.day == 1))

max_dec = np.max(FSM_dez.DATA).values
max_feb = np.max(FSM_feb.DATA).values
max_apr = np.max(FSM_apr.DATA).values

# make quick plots
fig, axs = plt.subplots(3, 4, figsize=[15, 7], frameon=False, subplot_kw={'projection': proj_map})
for tix in range(0, np.shape(axs)[1]):
    FSM_dez.isel(time=tix).DATA.plot(transform=proj_data, cmap='YlGnBu', add_colorbar=True, ax=axs[0, tix],
                                    cbar_kwargs={'label': "Snow Depth [m]"}, vmin=0, vmax=max_dec)
    FSM_feb.isel(time=tix).DATA.plot(transform=proj_data, cmap='YlGnBu', add_colorbar=True, ax=axs[1, tix],
                                    cbar_kwargs={'label': "Snow Depth [m]"}, vmin=0, vmax=max_feb)
    FSM_apr.isel(time=tix).DATA.plot(transform=proj_data, cmap='YlGnBu', add_colorbar=True, ax=axs[2, tix],
                                    cbar_kwargs={'label': "Snow Depth [m]"}, vmin=0, vmax=max_apr)
for ii in range(0, np.shape(axs)[0]):
    for jj in range(0, np.shape(axs)[1]):
        axs[ii, jj].axis('off')
        axs[ii, jj].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
plt.tight_layout()
plt.show()
fig.savefig(bf_out / 'comp_taylor' / Path('FSM' + '.png'), transparent=False, bbox_inches='tight')

FSM_feb_flat = FSM_feb.DATA.values.flatten()
FSM_feb_flat = FSM_feb_flat[~np.isnan(FSM_feb_flat)]
FSM_apr_flat = FSM_apr.DATA.values.flatten()
FSM_apr_flat = FSM_apr_flat[~np.isnan(FSM_apr_flat)]
FSM_dec_flat = FSM_dez.DATA.values.flatten()
FSM_dec_flat = FSM_dec_flat[~np.isnan(FSM_dec_flat)]

run_in_all_1km = ['CRUJRA_FILES_noLapse_OLD', 'CRUJRA_FILES_noLapse',
                  'CRUJRA_FILES_OLD', 'CRUJRA_FILES', 'OSHD_FILES_OLD', 'OSHD_FILES']

run_in_coarse = ['CRUJRA_FILES_05deg_cru_new', 'CRUJRA_FILES_05deg_cru_new_lapse', 'OSHD_FILES_05_new',
                 'CRUJRA_FILES_025deg_cru_new', 'CRUJRA_FILES_025deg_cru_new_lapse', 'OSHD_FILES_025_new']

label_in = ['FSM', 'Clim$_{CRU0.5^{\circ}}$+LU$_{Gl0.5^{\circ}}$', 'Clim$_{CRU^{*}0.5^{\circ}}$+LU$_{Gl0.5^{\circ}}$',
            'Clim$_{OSHD0.5^{\circ}}$+LU$_{Gl0.5^{\circ}}$',
            'Clim$_{CRU0.25^{\circ}}$+LU$_{Gl0.25^{\circ}}$', 'Clim$_{CRU^{*}0.25^{\circ}}$+LU$_{Gl0.25^{\circ}}$',
            'Clim$_{OSHD0.25^{\circ}}$+LU$_{Gl0.25^{\circ}}$',
            'Clim$_{CRU1km}$+LU$_{Gl1km}$', 'Clim$_{CRU1km}$+LU$_{HR1km}$',
            'Clim$_{CRU^{*}1km}$+LU$_{Gl1km}$', 'Clim$_{CRU^{*}1km}$+LU$_{HR1km}$',
            'Clim$_{OSHD1km}$+LU$_{Gl1km}$', 'Clim$_{OSHD1km}$+LU$_{HR1km}$'
            ]

run_in_all = run_in_coarse + run_in_all_1km

taylor_feb = []
taylor_dec = []
taylor_apr = []

sdev_feb = []
sdev_apr = []
sdev_dec = []
crmsd_feb = []
crmsd_apr = []
crmsd_dec = []
ccoef_feb = []
ccoef_apr = []
ccoef_dec = []

snow_in_025 = '*_025_SNOW_DEPTH*'

for run_in in run_in_all:
    if (run_in == 'CRUJRA_FILES_025deg_cru_new' or run_in == 'CRUJRA_FILES_025deg_cru_new_lapse' or
            run_in == 'OSHD_FILES_025_new'):
        snow_in_comp = xr.open_dataset(glob.glob(str(bf / 'OSHD_FILES' / snow_in_025))[0])
        snow_in = xr.open_dataset(glob.glob(str(bf / run_in / 'SNOW_DEPTH*'))[0])
        snow_in = snow_in.where(~np.isnan(snow_in_comp.DATA).values, np.nan)
    elif (run_in == 'CRUJRA_FILES_05deg_cru_new' or run_in == 'CRUJRA_FILES_05deg_cru_new_lapse' or
          run_in == 'OSHD_FILES_05_new'):
        snow_in_comp = xr.open_dataset(glob.glob(str(bf / 'OSHD_FILES' / snow_in_025))[0])
        snow_in = xr.open_dataset(glob.glob(str(bf / run_in / snow_in_025))[0])
        snow_in = snow_in.where(~np.isnan(snow_in_comp.DATA).values, np.nan)
    else:
        snow_in = xr.open_dataset(glob.glob(str(bf / run_in / snow_in_025))[0])
    clm5_feb = snow_in.isel(time=(snow_in.time.dt.month == 2) & (snow_in.time.dt.day == 1) &
                                 (snow_in.time.dt.year > 2015))
    clm5_apr = snow_in.isel(time=(snow_in.time.dt.month == 4) & (snow_in.time.dt.day == 1) &
                                 (snow_in.time.dt.year > 2015))
    clm5_dec = snow_in.isel(time=(snow_in.time.dt.month == 12) & (snow_in.time.dt.day == 1) &
                                 (snow_in.time.dt.year > 2014) & (snow_in.time.dt.year < 2019))

    max_dec = clm5_dec.DATA.max()
    max_feb = clm5_feb.DATA.max()
    max_apr = clm5_apr.DATA.max()

    #make quick plots
    fig, axs = plt.subplots(3, 4, figsize=[15, 7], frameon=False, subplot_kw={'projection': proj_map})
    for tix in range(0, np.shape(axs)[1]):
        clm5_dec.isel(time=tix).DATA.plot(transform=proj_data, cmap='YlGnBu', add_colorbar=True, ax=axs[0, tix],
                                        cbar_kwargs={'label': "Snow Depth [m]"}, vmin=0, vmax=max_dec)
        clm5_feb.isel(time=tix).DATA.plot(transform=proj_data, cmap='YlGnBu', add_colorbar=True, ax=axs[1, tix],
                                        cbar_kwargs={'label': "Snow Depth [m]"}, vmin=0, vmax=max_feb)
        clm5_apr.isel(time=tix).DATA.plot(transform=proj_data, cmap='YlGnBu', add_colorbar=True, ax=axs[2, tix],
                                        cbar_kwargs={'label': "Snow Depth [m]"}, vmin=0, vmax=max_apr)
    for ii in range(0, np.shape(axs)[0]):
        for jj in range(0, np.shape(axs)[1]):
            axs[ii, jj].axis('off')
            axs[ii, jj].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    plt.tight_layout()
    fig.savefig(bf_out / 'comp_taylor' / Path(run_in + '.png'), transparent=False, bbox_inches='tight')

    CLM5_feb_flat = clm5_feb.DATA.values.flatten()
    CLM5_feb_flat = CLM5_feb_flat[~np.isnan(CLM5_feb_flat)]

    CLM5_apr_flat = clm5_apr.DATA.values.flatten()
    CLM5_apr_flat = CLM5_apr_flat[~np.isnan(CLM5_apr_flat)]

    CLM5_dec_flat = clm5_dec.DATA.values.flatten()
    CLM5_dec_flat = CLM5_dec_flat[~np.isnan(CLM5_dec_flat)]

    stat_feb = sm.taylor_statistics(CLM5_feb_flat, FSM_feb_flat)
    stat_apr = sm.taylor_statistics(CLM5_apr_flat, FSM_apr_flat)
    stat_dec = sm.taylor_statistics(CLM5_dec_flat, FSM_dec_flat)

    taylor_feb.append(sm.taylor_statistics(CLM5_feb_flat, FSM_feb_flat))
    taylor_apr.append(sm.taylor_statistics(CLM5_apr_flat, FSM_apr_flat))
    taylor_dec.append(sm.taylor_statistics(CLM5_dec_flat, FSM_dec_flat))

    if run_in == 'CRUJRA_FILES_05deg_cru_new':
        sdev_feb.append(stat_feb['sdev'][0])
        crmsd_feb.append(stat_feb['crmsd'][0])
        ccoef_feb.append(stat_feb['ccoef'][0])

        sdev_apr.append(stat_apr['sdev'][0])
        crmsd_apr.append(stat_apr['crmsd'][0])
        ccoef_apr.append(stat_apr['ccoef'][0])

        sdev_dec.append(stat_dec['sdev'][0])
        crmsd_dec.append(stat_dec['crmsd'][0])
        ccoef_dec.append(stat_dec['ccoef'][0])

    sdev_feb.append(stat_feb['sdev'][1])
    crmsd_feb.append(stat_feb['crmsd'][1])
    ccoef_feb.append(stat_feb['ccoef'][1])

    sdev_apr.append(stat_apr['sdev'][1])
    crmsd_apr.append(stat_apr['crmsd'][1])
    ccoef_apr.append(stat_apr['ccoef'][1])

    sdev_dec.append(stat_dec['sdev'][1])
    crmsd_dec.append(stat_dec['crmsd'][1])
    ccoef_dec.append(stat_dec['ccoef'][1])

# Store statistics in arrays
sdev_feb = np.array(sdev_feb)
sdev_apr = np.array(sdev_apr)
sdev_dec = np.array(sdev_dec)

crmsd_feb = np.array(crmsd_feb)
crmsd_apr = np.array(crmsd_apr)
crmsd_dec = np.array(crmsd_dec)

ccoef_feb = np.array(ccoef_feb)
ccoef_apr = np.array(ccoef_apr)
ccoef_dec = np.array(ccoef_dec)


LEGEND_SUBPLOT = (0, 3)

SUBPLOTS_DATA = [
    {
        "axis_idx": (0, 0),
        "title": "(a) Early accumulation period",
        "y_label": True,
        "x_label": True,
        "observed": (sdev_dec[0], crmsd_dec[0], ccoef_dec[0]),
        "modeled": {
            label_in[1]: (sdev_dec[1], crmsd_dec[1], ccoef_dec[1]),
            label_in[2]: (sdev_dec[2], crmsd_dec[2], ccoef_dec[2]),
            label_in[3]: (sdev_dec[3], crmsd_dec[3], ccoef_dec[3]),
            label_in[4]: (sdev_dec[4], crmsd_dec[4], ccoef_dec[4]),
            label_in[5]: (sdev_dec[5], crmsd_dec[5], ccoef_dec[5]),
            label_in[6]: (sdev_dec[6], crmsd_dec[6], ccoef_dec[6]),
            label_in[7]: (sdev_dec[7], crmsd_dec[7], ccoef_dec[7]),
            label_in[8]: (sdev_dec[8], crmsd_dec[8], ccoef_dec[8]),
            label_in[9]: (sdev_dec[9], crmsd_dec[9], ccoef_dec[9]),
            label_in[10]: (sdev_dec[10], crmsd_dec[10], ccoef_dec[10]),
            label_in[11]: (sdev_dec[11], crmsd_dec[11], ccoef_dec[11]),
            label_in[12]: (sdev_dec[12], crmsd_dec[12], ccoef_dec[12])
        }
    }, {
        "axis_idx": (0, 1),
        "title": "(b) Mid-accumulation period",
        "y_label": False,
        "x_label": True,
        "observed": (sdev_feb[0], crmsd_feb[0], ccoef_feb[0]),
        "modeled": {
            label_in[1]: (sdev_feb[1], crmsd_feb[1], ccoef_feb[1]),
            label_in[2]: (sdev_feb[2], crmsd_feb[2], ccoef_feb[2]),
            label_in[3]: (sdev_feb[3], crmsd_feb[3], ccoef_feb[3]),
            label_in[4]: (sdev_feb[4], crmsd_feb[4], ccoef_feb[4]),
            label_in[5]: (sdev_feb[5], crmsd_feb[5], ccoef_feb[5]),
            label_in[6]: (sdev_feb[6], crmsd_feb[6], ccoef_feb[6]),
            label_in[7]: (sdev_feb[7], crmsd_feb[7], ccoef_feb[7]),
            label_in[8]: (sdev_feb[8], crmsd_feb[8], ccoef_feb[8]),
            label_in[9]: (sdev_feb[9], crmsd_feb[9], ccoef_feb[9]),
            label_in[10]: (sdev_feb[10], crmsd_feb[10], ccoef_feb[10]),
            label_in[11]: (ccoef_feb[11], ccoef_feb[11], ccoef_feb[11]),
            label_in[12]: (ccoef_feb[12], ccoef_feb[12], ccoef_feb[12])
        }
    }, {
        "axis_idx": (0, 2),
        "title": "(c) Ablation period",
        "y_label": False,
        "x_label": True,
        "observed": (sdev_apr[0], crmsd_apr[0], ccoef_apr[0]),
        "modeled": {
            label_in[1]: (sdev_apr[1], crmsd_apr[1], ccoef_apr[1]),
            label_in[2]: (sdev_apr[2], crmsd_apr[2], ccoef_apr[2]),
            label_in[3]: (sdev_apr[3], crmsd_apr[3], ccoef_apr[3]),
            label_in[4]: (sdev_apr[4], crmsd_apr[4], ccoef_apr[4]),
            label_in[5]: (sdev_apr[5], crmsd_apr[5], ccoef_apr[5]),
            label_in[6]: (sdev_apr[6], crmsd_apr[6], ccoef_apr[6]),
            label_in[7]: (sdev_apr[7], crmsd_apr[7], ccoef_apr[7]),
            label_in[8]: (sdev_apr[8], crmsd_apr[8], ccoef_apr[8]),
            label_in[9]: (sdev_apr[9], crmsd_apr[9], ccoef_apr[9]),
            label_in[10]: (sdev_apr[10], crmsd_apr[10], ccoef_apr[10]),
            label_in[11]: (ccoef_apr[11], ccoef_apr[11], ccoef_apr[11]),
            label_in[12]: (ccoef_apr[12], ccoef_apr[12], ccoef_apr[12])
        }
    }
]

existing_colors = [(217/255, 95/255, 2/255), (217/255, 95/255, 2/255), (117/255, 112/255, 179/255),
                   (1.0, 1.0, 1.0), (0.0, 0.0, 0.0)]
N = 6
col_in = distinctipy.get_colors(N, existing_colors, colorblind_type="Deuteranomaly")
# get colors which are as different as possible (also for colorblind)

MARKERS = {
    label_in[0]: {
        "marker": "^",
        "color_edge": "#000000",
        "color_face": "#000000",
        "markersize": 14
    },
    label_in[1]: {
        "marker": "o",
        "color_edge": 'grey',
        "color_face": col_in[0],
        "markersize": 9
    },
    label_in[2]: {
        "marker": "o",
        "color_edge": "grey",
        "color_face": col_in[1],
        "markersize": 9
    },
    label_in[3]: {
        "marker": "o",
        "color_edge": "grey",
        "color_face": col_in[2],
        "markersize": 9
    },
    label_in[4]: {
        "marker": "v",
        "color_edge": "grey",
        "color_face": col_in[3],
        "markersize": 9
    },
    label_in[5]: {
        "marker": "v",
        "color_edge": "grey",
        "color_face": col_in[4],
        "markersize": 9
    },
    label_in[6]: {
        "marker": "v",
        "color_edge": "grey",
        "color_face": col_in[5],
        "markersize": 9
    },
    label_in[7]: {
        "marker": "d",
        "color_edge": "grey",
        "color_face": (217/255, 95/255, 2/255),
        "markersize": 9
    },
    label_in[8]: {
        "marker": "*",
        "color_edge": "black",
        "color_face": (217/255, 95/255, 2/255),
        "markersize": 10
    },
    label_in[9]: {
        "marker": "d",
        "color_edge": "grey",
        "color_face": (117/255, 112/255, 179/255),
        "markersize": 9
    },
    label_in[10]: {
        "marker": "*",
        "color_edge": "black",
        "color_face": (117/255, 112/255, 179/255),
        "markersize": 10
    },
    label_in[11]: {
        "marker": "d",
        "color_edge": "grey",
        "color_face": (27/255, 158/255, 119/255),
        "markersize": 9
    },
    label_in[12]: {
        "marker": "*",
        "color_edge": "black",
        "color_face": (27/255, 158/255, 119/255),
        "markersize": 10
    }
}


# ## PLOT STYLE ################################################################# #

FONT_FAMILY = 'sans-serif'
FONT_SIZE = 10
# specify some styles for the correlation component
COLS_COR = {
    'grid': '#DDDDDD',
    'tick_labels': '#636161',
    'title': '#636161'
}

# specify some styles for the standard deviation
COLS_STD = {
    'grid': '#DDDDDD',
    'tick_labels': '#000000',
    'ticks': '#000000',
    'title': '#000000'
}

# specify some styles for the root mean square deviation
STYLES_RMS = {
    'color': '#9696e0',
    'linestyle': '--',
    'widthRMS': 1.8
}

plt.rcParams.update({'font.size': FONT_SIZE, 'font.family': FONT_FAMILY})

intervalsCOR = np.concatenate((np.arange(0, 1.0, 0.2), [0.9, 0.95, 0.99, 1]))


# create figure with 2 lines and 3 columns
fig_size = (12, 4)
fig, axs = plt.subplots(1, 4, figsize=fig_size)
del fig_size
# build subplot by subplot
for subplot_data in SUBPLOTS_DATA:

    # get subplot object and ensure it will be a square
    # y-axis labels will only appear on leftmost subplot
    print(subplot_data["axis_idx"])
    ax = axs[subplot_data["axis_idx"][1]]
    ax.set(adjustable='box', aspect='equal')
    amax = 1.1
    tick_rms = [0.2, 0.4, 0.6]
    tickSTD = [0.2, 0.4, 0.6, 0.8, 1]
    rincSTD = [0.2, 0.4, 0.6, 0.8, 1]
    title_ax_x = 0.04
    title_ax_y = -0.02
    if subplot_data["axis_idx"][1] == 0:
        amax = 0.5
        tick_rms = [0.2, 0.4]
        title_ax = 'FSM 1st Dec.'
        title_ax_x = 0.018
        title_ax_y = -0.01
        tickSTD = [0.2, 0.4]
        rincSTD = [0.2, 0.4]

    elif subplot_data["axis_idx"][1] == 1:
        title_ax = 'FSM 1st Feb.'
    else:
        title_ax = 'FSM 1st Apr.'

    # create the plot with the observed data
    stdev, crmsd, ccoef = subplot_data["observed"]
    sm.taylor_diagram(ax,
                      np.asarray((stdev, stdev)),
                      np.asarray((crmsd, crmsd)),
                      np.asarray((ccoef, ccoef)),
                      markercolors={
                          "face": MARKERS["FSM"]["color_edge"],
                          "edge": MARKERS["FSM"]["color_face"]
                      },
                      markersize=MARKERS["FSM"]["markersize"],
                      markersymbol=MARKERS["FSM"]["marker"],
                      styleOBS=':',
                      colOBS=MARKERS["FSM"]["color_edge"],
                      alpha=1.0,
                      titleSTD='off',
                      titleRMS='off',
                      showlabelsRMS='on',
                      tickRMS=tick_rms,
                      colRMS=STYLES_RMS['color'],
                      widthRMS=STYLES_RMS['widthRMS'],
                      tickRMSangle=115,
                      styleRMS=STYLES_RMS['linestyle'],
                      colscor=COLS_COR,
                      colsstd=COLS_STD,
                      styleCOR='-',
                      styleSTD=':',
                      tickCOR=intervalsCOR,
                      colframe='#DDDDDD',
                      labelweight='normal',
                      titlecorshape='linear',
                      axismax=amax,
                      tickSTD=tickSTD,
                      rincSTD=rincSTD)

    # add label below the marker
    ax.text(stdev+title_ax_x, title_ax_y, title_ax, verticalalignment="center",
            horizontalalignment="left",
            fontsize=FONT_SIZE + 1, fontweight="bold", rotation=0)

    # get rid of variables not to be used anymore
    del stdev, crmsd, ccoef

    # create one overlay for each model marker
    for model_id, (stdev, crmsd, ccoef) in subplot_data["modeled"].items():
        marker = MARKERS[model_id]
        sm.taylor_diagram(ax,
                          np.asarray((stdev, stdev)),
                          np.asarray((crmsd, crmsd)),
                          np.asarray((ccoef, ccoef)),
                          markersymbol=marker["marker"],
                          markercolors={
                              "face": marker["color_face"],
                              "edge": marker["color_edge"]
                          },
                          markersize=marker["markersize"],
                          alpha=1.0,
                          overlay='on',
                          styleCOR='-',
                          styleSTD='-')

        # get rid of variables not to be used anymore
        del model_id, stdev, crmsd, ccoef, marker

    # set titles (upper, left, bottom)
    ax.set_title(subplot_data["title"], loc="left", y=1.2, fontsize=FONT_SIZE + 3)

    # add y label
    if subplot_data["y_label"]:
        ax.set_ylabel("Standard Deviation [m]", fontfamily=FONT_FAMILY,
                      fontsize=FONT_SIZE + 3)

    # add xlabel or hide xticklabels
    if subplot_data["x_label"]:
        ax.set_xlabel("Standard Deviation [m]", fontfamily=FONT_FAMILY,
                      fontsize=FONT_SIZE + 3)
    else:
        ax.set_xticklabels(ax.get_xticklabels(), color=ax.get_facecolor())

    # just for the peace of mind...
    del subplot_data, ax

# create legend in the last subplot
ax = axs[LEGEND_SUBPLOT[1]]
ax.axis('off')

# build legend handles
legend_handles = []
legend_handles.append(mlines.Line2D([], [],
                                    color=STYLES_RMS['color'],
                                    linestyle=STYLES_RMS['linestyle'],
                                    label="RMSE"))

for marker_label, marker_desc in MARKERS.items():
    marker = mlines.Line2D([], [],
                           marker=marker_desc["marker"],
                           markersize=marker_desc["markersize"],
                           markerfacecolor=marker_desc["color_face"],
                           markeredgecolor=marker_desc["color_edge"],
                           linestyle='None',
                           label=marker_label)
    if marker_label != 'FSM':
        legend_handles.append(marker)
    del marker_label, marker_desc, marker

# create legend and free memory
ax.legend(handles=legend_handles, loc="center", fontsize=FONT_SIZE + 2)
del ax, legend_handles

# avoid some overlapping
plt.tight_layout()


bf_out = bf / 'figures_revisions'
fig.savefig(bf_out / 'taylor_upscaled.png', transparent=False, bbox_inches='tight')
