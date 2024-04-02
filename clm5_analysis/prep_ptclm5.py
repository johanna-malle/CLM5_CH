# -*- coding: utf-8 -*-
"""
Desc: quick script to convert PTCLM5 output to csv for further analysis
Created on 02.04.24 13:26
@author: malle
"""

import pandas as pd
from pathlib import Path
import xarray as xr
import glob
import os
import matplotlib.pyplot as plt


bf = Path('/home/lud11/malle/CLM5_CH/PTCLM5_nofor')
bf_clm5 = Path('/home/malle/CLM5_install/scratch/archive')

runs_all = ['PTCLM5_nofor_OSHD', 'PTCLM5_nofor_CRU_lapse', 'PTCLM5_nofor_CRU_Nolapse']

for run in runs_all:
    file_in_all = glob.glob(str(bf_clm5) + '/' + run + '/lnd/hist/*' + '.clm2.h0.*')
    path_out = bf / run
    path_out.mkdir(parents=True, exist_ok=True)
    for file_in in file_in_all:
        file_name = file_in.split('/')[-1].split('.')[0]
        file_name_out = path_out / Path(file_name + '_snow_related.csv')
        vars_in = xr.open_dataset(file_in)
        snow_depth = vars_in['SNOW_DEPTH'].data.squeeze()
        swe = vars_in['H2OSNO'].data.squeeze()
        rain = vars_in['RAIN'].data.squeeze()
        snow = vars_in['SNOW'].data.squeeze()
        fsno = vars_in['FSNO'].data.squeeze()
        time_clm5 = pd.to_datetime(vars_in.indexes['time'].to_datetimeindex())

        all_data = {'timeyears_clm': time_clm5, 'HS': snow_depth, 'SWE': swe, 'SCF': fsno, 'rain': rain, 'snow': snow}
        df = pd.DataFrame(all_data)
        df.set_index('timeyears_clm', inplace=True)
        df = df[df.index > '2014-01-01']
        df.to_csv(file_name_out)


path_FSM = '/home/lud11/malle/CLM5_CH/FSM_new/analysed_points'
all_files = glob.glob(os.path.join(path_FSM, "*.csv"))  # just do this once to get all ids
all_locs_comp = list((f.split('/')[-1]).split('_')[1] for f in all_files)
K = "MAE2"
K1 = "5DO"
all_locs = [i for i in all_locs_comp if (i != K and i != K1)]

path_oshd = bf / 'PTCLM5_nofor_OSHD'
path_cru_lapse = bf / 'PTCLM5_nofor_CRU_lapse'
path_cru_Nolapse = bf / 'PTCLM5_nofor_CRU_Nolapse'
bf_meas = Path('/home/lud11/malle/CLM5_CH/dvd_oshd')

bf_check = Path('/home/lud11/malle/CLM5_CH')
bf_oshd_luhr1km = bf_check / 'PTCLM_all_OSHD_newSurf'
bf_oshd_lugl1km = bf_check / 'PTCLM_all_OSHD_origSurf'

bf_plots = bf / 'plots'

for locs in all_locs:

    meas_in = pd.read_csv(glob.glob(os.path.join(bf_meas, "*" + locs + "*.csv"))[0]).set_index('time_HS_meas')
    meas_in = meas_in.loc[~meas_in.index.duplicated(keep='first')]  # this is necessary since some seasons overlapped..

    new_cru_lapse = pd.read_csv(glob.glob(os.path.join(path_cru_lapse,
                                                          "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')
    new_cru_nolapse = pd.read_csv(glob.glob(os.path.join(path_cru_Nolapse,
                                                          "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')
    new_oshd = pd.read_csv(glob.glob(os.path.join(path_oshd,
                                                          "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')

    old_oshd_hr = pd.read_csv(glob.glob(os.path.join(bf_oshd_luhr1km,
                                                          "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')
    old_oshd_gl = pd.read_csv(glob.glob(os.path.join(bf_oshd_lugl1km,
                                                          "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')

    fig = plt.figure(figsize=(12, 13))
    axes = fig.add_subplot(311)
    axes.yaxis.grid(True)
    axes.plot(pd.to_datetime(new_cru_nolapse.index), new_cru_nolapse.HS, label=r'Clim$_{CRU}$')
    axes.plot(pd.to_datetime(new_cru_lapse.index), new_cru_lapse.HS, label=r'Clim$_{CRU*}$')
    axes.plot(pd.to_datetime(new_oshd.index), new_oshd.HS, label=r'Clim$_{OSHD}$')
    axes.plot(pd.to_datetime(old_oshd_gl.index), old_oshd_gl.HS, label=r'Clim$_{OSHD}$LU$_{Gl1KM}$')
    axes.plot(pd.to_datetime(old_oshd_hr.index), old_oshd_hr.HS, label=r'Clim$_{OSHD}$LU$_{HR1KM}$')
    axes.plot(pd.to_datetime(meas_in.index), meas_in.HS_meas, color='k', linestyle='--', linewidth=2, label='MEASURED')
    axes.set_xlim([pd.to_datetime('2015-09-01'), pd.to_datetime('2020-07-01')])
    plt.legend()
    axes.set_ylabel('SNOW DEPTH [m]')

    axes = fig.add_subplot(312)
    axes.yaxis.grid(True)
    axes.plot(pd.to_datetime(new_cru_nolapse.index), new_cru_nolapse.SCF, label=r'Clim$_{CRU}$')
    axes.plot(pd.to_datetime(new_cru_lapse.index), new_cru_lapse.SCF, label=r'Clim$_{CRU*}$')
    axes.plot(pd.to_datetime(new_oshd.index), new_oshd.SCF, label=r'Clim$_{OSHD}$')
    axes.set_xlim([pd.to_datetime('2015-09-01'), pd.to_datetime('2020-07-01')])
    axes.set_ylabel('SCF [-]')

    axes = fig.add_subplot(313)
    axes.yaxis.grid(True)
    axes.plot(pd.to_datetime(new_cru_nolapse.index), new_cru_nolapse.snow, label=r'Clim$_{CRU}$')
    axes.plot(pd.to_datetime(new_cru_lapse.index), new_cru_lapse.snow, label=r'Clim$_{CRU*}$')
    axes.plot(pd.to_datetime(new_oshd.index), new_oshd.snow, label=r'Clim$_{OSHD}$')
    axes.set_xlim([pd.to_datetime('2015-09-01'), pd.to_datetime('2020-07-01')])
    axes.set_ylabel('atmospheric snow [mm/s]')

    plt.tight_layout()
    file_out = bf_plots / Path('snow_comp_'+locs+'.png')
    fig.savefig(file_out, facecolor='white', transparent=False, bbox_inches='tight')








