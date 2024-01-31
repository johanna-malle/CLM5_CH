# -*- coding: utf-8 -*-
"""
Desc: 
Created on 30.01.24 10:27
@author: malle
"""

import pandas as pd
from pathlib import Path
import xarray as xr
import glob
import argparse
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import skill_metrics as sm
import numpy as np

# make switch to windows if working from home
mach = 'windows'
if mach == 'linux':
    bf = Path('/home/lud11/malle/CLM5_CH')
else:
    bf = Path('L:\malle\CLM5_CH')

bf_oshd = bf / 'FSM_new/analysed_grid_python_out'
all_FSM = xr.open_mfdataset(glob.glob(str(bf_oshd / '*025*')))


FSM_feb = all_FSM.isel(time=(all_FSM.time.dt.month == 2) & (all_FSM.time.dt.day == 1))
FSM_apr = all_FSM.isel(time=(all_FSM.time.dt.month == 4) & (all_FSM.time.dt.day == 1))
FSM_dez = all_FSM.isel(time=(all_FSM.time.dt.month == 12) & (all_FSM.time.dt.day == 1))

FSM_feb_flat = FSM_feb.DATA.values.flatten()
FSM_feb_flat = FSM_feb_flat[~np.isnan(FSM_feb_flat)]
FSM_apr_flat = FSM_apr.DATA.values.flatten()
FSM_apr_flat = FSM_apr_flat[~np.isnan(FSM_apr_flat)]
FSM_dec_flat = FSM_dez.DATA.values.flatten()
FSM_dec_flat = FSM_dec_flat[~np.isnan(FSM_dec_flat)]

run_in_all_1km = ['OSHD_FILES', 'OSHD_FILES_OLD', 'CRUJRA_FILES_noLapse', 'CRUJRA_FILES_noLapse_OLD', 'CRUJRA_FILES',
                  'CRUJRA_FILES_OLD']

run_in_coarse = ['CRUJRA_FILES_025deg_cru_new', 'CRUJRA_FILES_05deg_cru_new',
                 'CRUJRA_FILES_025deg_cru_new_lapse', 'CRUJRA_FILES_05deg_cru_new_lapse']

run_in_all = run_in_all_1km + run_in_coarse

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

for run_in in run_in_all:
    if run_in == 'CRUJRA_FILES_025deg_cru_new' or  run_in =='CRUJRA_FILES_025deg_cru_new_lapse' :
        snow_in_comp = xr.open_dataset(glob.glob(str(bf / 'OSHD_FILES' / '*_025_SNOW_DEPTH*'))[0])
        snow_in = xr.open_dataset(glob.glob(str(bf / run_in / 'SNOW_DEPTH*'))[0])
        snow_in = snow_in.where(~np.isnan(snow_in_comp.DATA).values, np.nan)
    elif run_in == 'CRUJRA_FILES_05deg_cru_new' or  run_in =='CRUJRA_FILES_05deg_cru_new_lapse' :
        snow_in_comp = xr.open_dataset(glob.glob(str(bf / 'OSHD_FILES' / '*_025_SNOW_DEPTH*'))[0])
        snow_in = xr.open_dataset(glob.glob(str(bf / run_in / '*_025_SNOW_DEPTH*'))[0])
        snow_in = snow_in.where(~np.isnan(snow_in_comp.DATA).values, np.nan)
    else:
        snow_in = xr.open_dataset(glob.glob(str(bf / run_in / '*_025_SNOW_DEPTH*'))[0])
    clm5_feb = snow_in.isel(time=(snow_in.time.dt.month == 2) & (snow_in.time.dt.day == 1) & (snow_in.time.dt.year > 2015))
    clm5_apr = snow_in.isel(time=(snow_in.time.dt.month == 4) & (snow_in.time.dt.day == 1) & (snow_in.time.dt.year > 2015))
    clm5_dec = snow_in.isel(time=(snow_in.time.dt.month == 12) & (snow_in.time.dt.day == 1)
                                 & (snow_in.time.dt.year > 2014) & (snow_in.time.dt.year < 2019))

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

    if run_in == 'OSHD_FILES':
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


LEGEND_SUBPLOT = (1, 2)

SUBPLOTS_DATA = [
    {
        "axis_idx": (0, 0),
        "title": "(a) 1 hour lead time",
        "y_label": True,
        "x_label": False,
        "observed": (29.91, 0.00, 1.00),
        "modeled": {
            "Model \u03B1": (30.28, 16.84, 0.84),
            "Model \u03B2": (28.54,  8.83, 0.96),
            "Model \u03B4": (28.02,  9.02, 0.95),
            "Model \u03B8": (28.51,  8.44, 0.96),
            "Model \u03C1": (28.39,  8.47, 0.96)
        }
    }, {
        "axis_idx": (0, 1),
        "title": "(b) 2 hours lead time",
        "y_label": False,
        "x_label": False,
        "observed": (29.91, 0.00, 1.00),
        "modeled": {
            "Model \u03B1": (25.45, 17.58, 0.81),
            "Model \u03B2": (26.07, 14.69, 0.87),
            "Model \u03B4": (24.57, 16.08, 0.84),
            "Model \u03B8": (26.43, 15.26, 0.86),
            "Model \u03C1": (26.57, 13.61, 0.89)
        }
    }, {
        "axis_idx": (0, 2),
        "title": "(c) 3 hours lead time",
        "y_label": False,
        "x_label": True,
        "observed": (29.91, 0.00, 1.00),
        "modeled": {
            "Model \u03B1": (19.37, 23.33, 0.64),
            "Model \u03B2": (21.98, 21.00, 0.72),
            "Model \u03B4": (20.02, 21.74, 0.70),
            "Model \u03B8": (22.46, 20.37, 0.74),
            "Model \u03C1": (23.32, 19.30, 0.77)
        }
    }, {
        "axis_idx": (1, 0),
        "title": "(d) 4 hours lead time",
        "y_label": True,
        "x_label": True,
        "observed": (29.91, 0.00, 1.00),
        "modeled": {
            "Model \u03B1": (15.46, 27.45, 0.46),
            "Model \u03B2": (16.64, 25.49, 0.56),
            "Model \u03B4": (14.74, 25.94, 0.54),
            "Model \u03B8": (19.71, 25.73, 0.56),
            "Model \u03C1": (19.30, 23.28, 0.66)
        }
    }, {
        "axis_idx": (1, 1),
        "title": "(e) 5 hours lead time",
        "y_label": True,
        "x_label": True,
        "observed": (29.91, 0.00, 1.00),
        "modeled": {
            "Model \u03B1": (10.33, 30.00, 0.28),
            "Model \u03B2": (16.54, 30.09, 0.33),
            "Model \u03B4": (12.03, 28.26, 0.43),
            "Model \u03B8": (15.76, 29.92, 0.34),
            "Model \u03C1": (15.87, 27.06, 0.50)
        }
    }
]

MARKERS = {
    "Observed": {
        "marker": "^",
        "color_edge": "#000000",
        "color_face": "#000000",
        "markersize": 9
    },
    "Model \u03B1": {
        "marker": "o",
        "color_edge": "#000000",
        "color_face": "#777777",
        "markersize": 9
    },
    "Model \u03B2": {
        "marker": "D",
        "color_edge": "#AA0000",
        "color_face": "#DD3333",
        "markersize": 9
    },
    "Model \u03B4": {
        "marker": "v",
        "color_edge": "#00AA00",
        "color_face": "#33DD33",
        "markersize": 9
    },
    "Model \u03B8": {
        "marker": "s",
        "color_edge": "#0000AA",
        "color_face": "#3333DD",
        "markersize": 9
    },
    "Model \u03C1": {
        "marker": "*",
        "color_edge": "#D4AF37",
        "color_face": "#FFD700",
        "markersize": 12
    }
}


# ## PLOT STYLE ################################################################# #

FONT_FAMILY = 'Times New Roman'
FONT_SIZE = 9

label_in = ['FSM','OSHD_FILES', 'OSHD_FILES_OLD', 'CRUJRA_FILES_noLapse', 'CRUJRA_FILES_noLapse_OLD', 'CRUJRA_FILES',
                  'CRUJRA_FILES_OLD']

plt.figure(num=1, figsize=(8, 6))
sm.taylor_diagram(sdev_feb,crmsd_feb,ccoef_feb,
                  markerLabel=label_in, markerLabelColor='r', markerLegend='on', markerColor ='r', markerSize = 6,
                  styleOBS='-', colOBS='g', markerobs='o', colRMS='m', styleRMS=':', widthRMS=2.0,
                  tickRMSangle=115, showlabelsRMS='on', tickRMS = [0.2, 0.4, 0.6], colCOR='k', styleCOR='--', widthCOR = 1.0,
                  titleOBS='FSM2 1st Feb.')

plt.show()

run_in_05 = ['OSHD_FILES_05_new', 'CRUJRA_FILES_05deg_cru_new', 'CRUJRA_FILES_05deg_cru_new_lapse']
run_in_025 = ['OSHD_FILES_025_new', 'CRUJRA_FILES_025deg_cru_new', 'CRUJRA_FILES_025deg_cru_new_lapse']


