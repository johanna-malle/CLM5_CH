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
import xesmf as xe
import scipy.io
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from skill_metrics import centered_rms_dev
import skill_metrics as sm


# make switch to windows if working from home
mach = 'linux'
if mach == 'linux':
    bf = Path('/home/lud11/malle/CLM5_CH')
else:
    bf = Path('L:\malle\CLM5_CH')

bf_oshd = bf / 'FSM_new/analysed_grid_python_out'
all_FSM = xr.open_mfdataset(glob.glob(str(bf_oshd) + '/*025*'))
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

for run_in in run_in_all_1km:
    snow_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/upscale_025_SNOW_DEPTH*')[0])
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




label_in = ['FSM','OSHD_FILES', 'OSHD_FILES_OLD', 'CRUJRA_FILES_noLapse', 'CRUJRA_FILES_noLapse_OLD', 'CRUJRA_FILES',
                  'CRUJRA_FILES_OLD']

plt.figure(num=1, figsize=(8, 6))
sm.taylor_diagram(sdev_feb,crmsd_feb,ccoef_feb,
                  markerLabel = label_in, markerLabelColor = 'r', markerLegend = 'on', markerColor = 'r', markerSize = 6,
                  styleOBS = '-', colOBS = 'g', markerobs = 'o', colRMS = 'm', styleRMS = ':', widthRMS = 2.0,
                  tickRMSangle = 115, showlabelsRMS = 'on', tickRMS = [0.2, 0.4, 0.6],colCOR = 'k', styleCOR = '--', widthCOR = 1.0,
                  titleOBS = 'FSM2 1st Feb.')

plt.show()

run_in_05 = ['OSHD_FILES_05_new', 'CRUJRA_FILES_05deg_cru_new', 'CRUJRA_FILES_05deg_cru_new_lapse']
run_in_025 = ['OSHD_FILES_025_new', 'CRUJRA_FILES_025deg_cru_new', 'CRUJRA_FILES_025deg_cru_new_lapse']


