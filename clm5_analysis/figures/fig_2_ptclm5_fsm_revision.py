# -*- coding: utf-8 -*-
"""
Desc: Script to produce Figure 2.
"Comparisons of point-scale model simulations to observations of snow depth (HS)
across all simulated snow seasons (October–July)"
Created on 08.12.22 18:02
@author: malle
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import os
import seaborn as sns
import matplotlib
import glob
import platform


def mae(y_true, predictions):
    y_true, predictions = np.array(y_true), np.array(predictions)
    return np.round(np.mean(np.abs(y_true - predictions)), 2)


def rmse(y_true, predictions):
    y_true, predictions = np.array(y_true), np.array(predictions)
    return np.round(np.sqrt(((predictions - y_true) ** 2).mean()), 2)


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

if platform.system() == 'Windows':
    bf = Path('L:\malle\CLM5_CH')
else:
    bf = Path('/home/lud11/malle/CLM5_CH')

path_FSM = bf / 'FSM_new' / 'analysed_points'
all_files = glob.glob(os.path.join(path_FSM, "*.csv"))  # just do this once to get all ids

if platform.system() == 'Windows':
    all_locs_comp = list((f.split('\\')[-1]).split('_')[1] for f in all_files)
else:
    all_locs_comp = list((f.split('/')[-1]).split('_')[1] for f in all_files)

K = "MAE2"
K1 = "5DO"
all_locs = [i for i in all_locs_comp if (i != K and i != K1)]
path_oshd_newSurf = bf / 'PTCLM_all_OSHD_newSurf'
path_oshd_origSurf = bf / 'PTCLM_all_OSHD_origSurf'
path_crujra_newSurf = bf / 'PTCLM_all_CRUJRA_Nolapse_newSurf'
path_crujra_origSurf = bf / 'PTCLM_all_CRUJRA_Nolapse_origSurf'
path_crujraP_newSurf = bf / 'PTCLM_all_CRUJRA_lapse_newSurf'
path_crujraP_origSurf = bf / 'PTCLM_all_CRUJRA_lapse_origSurf'
path_crujraP_nofor = bf / 'PTCLM5_nofor' / 'PTCLM5_nofor_CRU_lapse'
path_crujra_nofor = bf / 'PTCLM5_nofor' / 'PTCLM5_nofor_CRU_Nolapse'
path_oshd_nofor = bf / 'PTCLM5_nofor' / 'PTCLM5_nofor_OSHD'

bf_meas = bf / 'dvd_oshd'
all_data_1000 = pd.DataFrame([])
all_data_2000 = pd.DataFrame([])
all_data_3000 = pd.DataFrame([])

# Read each CSV file in dir "path/to/root_dir"
all_elev = glob.glob(os.path.join(path_FSM, "*.csv"))
if platform.system() == 'Windows':
    all_names_elev = list((f.split('\\')[-1]).split('_')[1] for f in all_elev)
    all_elev_elev = list((f.split('\\')[-1]).split('_')[-1].split('.')[0] for f in all_elev)
else:
    all_names_elev = list((f.split('/')[-1]).split('_')[1] for f in all_elev)
    all_elev_elev = list((f.split('/')[-1]).split('_')[-1].split('.')[0] for f in all_elev)
all_elev_elev = list(map(int, all_elev_elev))
elev_comp = list(zip(all_names_elev, all_elev_elev))

# analysis per elevation band
less_1000 = [x[0] for x in elev_comp if x[1] <= 1000]
bw_1000_2000 = [x[0] for x in elev_comp if (x[1] > 1000 and x[1] < 2000)]
above_2000 = [x[0] for x in elev_comp if x[1] >= 2000]

for locs in all_locs:
    meas_in = pd.read_csv(glob.glob(os.path.join(bf_meas, "*" + locs + "*.csv"))[0]).set_index('time_HS_meas')
    meas_in = meas_in.loc[~meas_in.index.duplicated(keep='first')]  # this is necessary since some seasons overlapped..
    meas_in.set_index(pd.to_datetime(meas_in.index), inplace=True)
    meas_in.dropna(inplace=True)

    print(locs)
    name_out = bf / 'snow_comp' / Path('station_'+locs+'_nov_may.csv')
    name_out_diff = bf / 'snow_comp' / Path('station_'+locs+'_diff_nov_may.csv')

    oshd_newSurf = pd.read_csv(glob.glob(os.path.join(path_oshd_newSurf,
                                                      "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')
    oshd_origSurf = pd.read_csv(glob.glob(os.path.join(path_oshd_origSurf,
                                                       "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')
    crujra_newSurf = pd.read_csv(glob.glob(os.path.join(path_crujra_newSurf,
                                                        "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')
    crujra_origSurf = pd.read_csv(glob.glob(os.path.join(path_crujra_origSurf,
                                                         "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')
    crujraP_newSurf = pd.read_csv(glob.glob(os.path.join(path_crujraP_newSurf,
                                                         "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')
    crujraP_origSurf = pd.read_csv(glob.glob(os.path.join(path_crujraP_origSurf,
                                                          "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')
    crujraP_nofor = pd.read_csv(glob.glob(os.path.join(path_crujraP_nofor,
                                                       "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')
    crujra_nofor = pd.read_csv(glob.glob(os.path.join(path_crujra_nofor,
                                                      "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')
    oshd_nofor = pd.read_csv(glob.glob(os.path.join(path_oshd_nofor,
                                                    "*"+locs+"*.csv"))[0]).set_index('timeyears_clm')

    oshd_newSurf.set_index(pd.to_datetime(oshd_newSurf.index), inplace=True)
    oshd_origSurf.set_index(pd.to_datetime(oshd_origSurf.index), inplace=True)
    oshd_nofor.set_index(pd.to_datetime(oshd_nofor.index), inplace=True)

    crujra_newSurf.set_index(pd.to_datetime(crujra_newSurf.index), inplace=True)
    crujra_origSurf.set_index(pd.to_datetime(crujra_origSurf.index), inplace=True)
    crujra_nofor.set_index(pd.to_datetime(crujra_nofor.index), inplace=True)

    crujraP_newSurf.set_index(pd.to_datetime(crujraP_newSurf.index), inplace=True)
    crujraP_origSurf.set_index(pd.to_datetime(crujraP_origSurf.index), inplace=True)
    crujraP_nofor.set_index(pd.to_datetime(crujraP_nofor.index), inplace=True)

    result = pd.concat([oshd_newSurf[['HS', 'SWE']], oshd_origSurf[['HS', 'SWE']], oshd_nofor[['HS', 'SWE']],
                        crujra_newSurf[['HS', 'SWE']], crujra_origSurf[['HS', 'SWE']], crujra_nofor[['HS', 'SWE']],
                        crujraP_newSurf[['HS', 'SWE']], crujraP_origSurf[['HS', 'SWE']],
                        crujraP_nofor[['HS', 'SWE']]], axis=1, join='inner',
                       keys=('oshd_newSurf', 'oshd_origSurf', 'oshd_nofor',
                             'crujra_newSurf', 'crujra_origSurf', 'crujra_nofor',
                             'crujraP_newSurf', 'crujraP_origSurf', 'crujraP_nofor'))

    # just do this once to get all ids
    FSM = pd.read_csv(glob.glob(os.path.join(path_FSM, "*"+locs+"*.csv"))[0])  # .set_index('time_stamps_jim')
    FSM.loc[:, 'time_stamp_comp'] = FSM.time_stamps_jim.map(lambda x: x[:-9])
    FSM.set_index('time_stamp_comp', inplace=True)
    FSM.drop(['time_stamps_jim', 'scf_jim', 'SWE_jim'], axis=1, inplace=True)
    time_filter = pd.to_datetime(FSM.index)
    FSM1 = FSM[((time_filter.month < 8) | (time_filter.month > 9)) & (time_filter.year < 2021)]  # only comp. oct-july
    FSM1.set_index(pd.to_datetime(FSM1.index), inplace=True)

    result_all = pd.concat([result, FSM1], axis=1, join='inner')
    result_all.to_csv(name_out)
    result_all_meas = pd.concat([result_all, meas_in], axis=1, join='inner')

    if locs in less_1000:
        all_data_1000 = pd.concat([all_data_1000, result_all_meas])
    elif locs in bw_1000_2000:
        all_data_2000 = pd.concat([all_data_2000, result_all_meas])
    else:
        all_data_3000 = pd.concat([all_data_3000, result_all_meas])

    diff_all = result_all_meas.copy()
    diff_all['diff_oshd_new'] = diff_all[('oshd_newSurf', 'HS')].sub(diff_all['HS_meas'], axis=0)
    diff_all['diff_oshd_orig'] = diff_all[('oshd_origSurf', 'HS')].sub(diff_all['HS_meas'], axis=0)
    diff_all['diff_oshd_nofor'] = diff_all[('oshd_nofor', 'HS')].sub(diff_all['HS_meas'], axis=0)
    diff_all['diff_cru_new'] = diff_all[('crujra_newSurf', 'HS')].sub(diff_all['HS_meas'], axis=0)
    diff_all['diff_cru_orig'] = diff_all[('crujra_origSurf', 'HS')].sub(diff_all['HS_meas'], axis=0)
    diff_all['diff_cru_nofor'] = diff_all[('crujra_nofor', 'HS')].sub(diff_all['HS_meas'], axis=0)
    diff_all['diff_cruP_new'] = diff_all[('crujraP_newSurf', 'HS')].sub(diff_all['HS_meas'], axis=0)
    diff_all['diff_cruP_orig'] = diff_all[('crujraP_origSurf', 'HS')].sub(diff_all['HS_meas'], axis=0)
    diff_all['diff_cruP_nofor'] = diff_all[('crujraP_nofor', 'HS')].sub(diff_all['HS_meas'], axis=0)
    diff_all['FSM'] = diff_all['HS_jim'].sub(diff_all['HS_meas'], axis=0)

    diff_all.drop(columns=[('oshd_newSurf', 'HS'), ('oshd_newSurf', 'SWE'),
                           ('oshd_origSurf', 'HS'), ('oshd_origSurf', 'SWE'),
                           ('oshd_nofor', 'HS'), ('oshd_nofor', 'SWE'),
                           ('crujra_newSurf', 'HS'), ('crujra_newSurf', 'SWE'),
                           ('crujra_origSurf', 'HS'), ('crujra_origSurf', 'SWE'),
                           ('crujra_nofor', 'HS'), ('crujra_nofor', 'SWE'),
                           ('crujraP_newSurf', 'HS'), ('crujraP_newSurf', 'SWE'),
                           ('crujraP_origSurf', 'HS'), ('crujraP_origSurf', 'SWE'),
                           ('crujraP_nofor', 'HS'), ('crujraP_nofor', 'SWE'), 'HS_jim'], inplace=True)

    diff_all.to_csv(name_out_diff)

all_data_1000.drop(columns=[('oshd_newSurf', 'SWE'), ('oshd_origSurf', 'SWE'), ('oshd_nofor', 'SWE'),
                   ('crujra_newSurf', 'SWE'), ('crujra_origSurf', 'SWE'), ('crujra_nofor', 'SWE'),
                            ('crujraP_newSurf', 'SWE'), ('crujraP_origSurf', 'SWE'), ('crujraP_nofor', 'SWE')],
                   inplace=True)

all_data_2000.drop(columns=[('oshd_newSurf', 'SWE'), ('oshd_origSurf', 'SWE'), ('oshd_nofor', 'SWE'),
                   ('crujra_newSurf', 'SWE'), ('crujra_origSurf', 'SWE'), ('crujra_nofor', 'SWE'),
                            ('crujraP_newSurf', 'SWE'), ('crujraP_origSurf', 'SWE'), ('crujraP_nofor', 'SWE')],
                   inplace=True)

all_data_3000.drop(columns=[('oshd_newSurf', 'SWE'), ('oshd_origSurf', 'SWE'), ('oshd_nofor', 'SWE'),
                   ('crujra_newSurf', 'SWE'), ('crujra_origSurf', 'SWE'), ('crujra_nofor', 'SWE'),
                            ('crujraP_newSurf', 'SWE'), ('crujraP_origSurf', 'SWE'), ('crujraP_nofor', 'SWE')],
                   inplace=True)
all_data_1000_nona = all_data_1000.dropna()
all_data_2000_nona = all_data_2000.dropna()
all_data_3000_nona = all_data_3000.dropna()

rmse_1000 = []
mae_1000 = []
rmse_2000 = []
mae_2000 = []
rmse_3000 = []
mae_3000 = []
col_names = []

for col_name in all_data_1000_nona.keys():
    rmse_1000.append(rmse(all_data_1000_nona.HS_meas, all_data_1000_nona[col_name]))
    mae_1000.append(mae(all_data_1000_nona.HS_meas, all_data_1000_nona[col_name]))
    rmse_2000.append(rmse(all_data_2000_nona.HS_meas, all_data_2000_nona[col_name]))
    mae_2000.append(mae(all_data_2000_nona.HS_meas, all_data_2000_nona[col_name]))
    rmse_3000.append(rmse(all_data_3000_nona.HS_meas, all_data_3000_nona[col_name]))
    mae_3000.append(mae(all_data_3000_nona.HS_meas, all_data_3000_nona[col_name]))
    if np.size(col_name) == 1:
        col_names.append(col_name)
    else:
        col_names.append(col_name[0])

rmse_1000 = pd.DataFrame(rmse_1000, index=col_names).T
rmse_2000 = pd.DataFrame(rmse_2000, index=col_names).T
rmse_3000 = pd.DataFrame(rmse_3000, index=col_names).T

mae_1000 = pd.DataFrame(mae_1000, index=col_names).T
mae_2000 = pd.DataFrame(mae_2000, index=col_names).T
mae_3000 = pd.DataFrame(mae_3000, index=col_names).T

dfs = []
dfs_1000 = []
dfs_2000 = []
dfs_3000 = []
for file in Path(bf / 'snow_comp').glob("*diff_nov_may.csv"):
    dfs.append(pd.read_csv(file))

    if any(ext in str(file) for ext in less_1000):
        dfs_1000.append(pd.read_csv(file))
    elif any(ext in str(file) for ext in bw_1000_2000):
        dfs_2000.append(pd.read_csv(file))
    elif any(ext in str(file) for ext in above_2000):
        dfs_3000.append(pd.read_csv(file))

df = pd.concat(dfs)  # put the dataframes to a single dataframe
df_1000 = pd.concat(dfs_1000)
df_2000 = pd.concat(dfs_2000)
df_3000 = pd.concat(dfs_3000)

plt.rc('font', family='serif')

fig = plt.figure(figsize=(12, 12))
my_pal = {"diff_cru_orig": (217/255, 95/255, 2/255), "diff_cru_new": (217/255, 95/255, 2/255),
          "diff_cru_nofor": (217/255, 95/255, 2/255),
          "diff_cruP_orig": (117/255, 112/255, 179/255), "diff_cruP_new": (117/255, 112/255, 179/255),
          "diff_cruP_nofor": (117/255, 112/255, 179/255),
          "diff_oshd_orig": (27/255, 158/255, 119/255), "diff_oshd_new": (27/255, 158/255, 119/255),
          "diff_oshd_nofor": (27/255, 158/255, 119/255), "FSM": "dimgray"}

flierprops = dict(marker='o', markersize=4.5, markeredgecolor='gray', markerfacecolor='silver', alpha=0.45)
axes = fig.add_subplot(311)
ax = sns.boxplot(data=(df_1000[['diff_cru_orig', 'diff_cru_new', 'diff_cru_nofor', 'diff_cruP_orig', 'diff_cruP_new',
                                'diff_cruP_nofor', 'diff_oshd_orig', 'diff_oshd_new', 'diff_oshd_nofor', 'FSM']]),
                 ax=axes, palette=my_pal, flierprops=flierprops, linewidth=1.6, saturation=0.8, showfliers=False,
                 meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": "6"},
                 medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), showmeans=True)
axes.set_title(r"$\bf{(a)}$"+' Locations < 1000m (n = '+str(len(less_1000))+')', loc='left', fontsize=16)
axes.set_ylim([-2.1, 2.1])
axes.yaxis.grid(True)
axes.set_xticklabels([])
axes.set_ylabel(u'Δ HS [m] ', fontsize=14)
ax.axhline(0, color='k', linewidth=2.5)

for label in (axes.get_xticklabels() + axes.get_yticklabels()):
    label.set_fontsize(13)

for patch in ax.patches[1:-1:3]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.6))

for patch in ax.patches[2:-1:3]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.2))

plt.text(0.015, 0.08, 'RMSE='+str(rmse_1000.crujra_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.015, 0.02, 'MAE='+str(mae_1000.crujra_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.115, 0.08, 'RMSE='+str(rmse_1000.crujra_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.115, 0.02, 'MAE='+str(mae_1000.crujra_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.22, 0.08, 'RMSE='+str(rmse_1000.crujra_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.22, 0.02, 'MAE='+str(mae_1000.crujra_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.32, 0.08, 'RMSE='+str(rmse_1000.crujraP_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.32, 0.02, 'MAE='+str(mae_1000.crujraP_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.42, 0.08, 'RMSE='+str(rmse_1000.crujraP_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.42, 0.02, 'MAE='+str(mae_1000.crujraP_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.52, 0.08, 'RMSE='+str(rmse_1000.crujraP_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.52, 0.02, 'MAE='+str(mae_1000.crujraP_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.625, 0.08, 'RMSE='+str(rmse_1000.oshd_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.625, 0.02, 'MAE='+str(mae_1000.oshd_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.725, 0.08, 'RMSE='+str(rmse_1000.oshd_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.725, 0.02, 'MAE='+str(mae_1000.oshd_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.82, 0.08, 'RMSE='+str(rmse_1000.oshd_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.82, 0.02, 'MAE='+str(mae_1000.oshd_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.915, 0.08, 'RMSE='+str(rmse_1000.HS_jim.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.915, 0.02, 'MAE='+str(mae_1000.HS_jim.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

axes = fig.add_subplot(312)
ax = sns.boxplot(data=(df_2000[['diff_cru_orig', 'diff_cru_new', 'diff_cru_nofor', 'diff_cruP_orig', 'diff_cruP_new',
                                'diff_cruP_nofor', 'diff_oshd_orig', 'diff_oshd_new', 'diff_oshd_nofor', 'FSM']]),
                 ax=axes, palette=my_pal, flierprops=flierprops, linewidth=1.6, saturation=0.8, showfliers=False,
                 meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": "6"},
                 medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), showmeans=True)
ax.axhline(0, color='k', linewidth=2.5)
axes.set_ylim([-2.1, 2.1])

axes.set_title(r"$\bf{(b)}$"+' Locations 1000-2000m (n = '+str(len(bw_1000_2000))+')', loc='left', fontsize=16)
axes.yaxis.grid(True)
axes.set_ylabel(u'Δ HS [m] ', fontsize=14)
axes.set_xticklabels([])
for label in (axes.get_xticklabels() + axes.get_yticklabels()):
    label.set_fontsize(13)

for patch in ax.patches[1:-1:3]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.6))

for patch in ax.patches[2:-1:3]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.2))

plt.text(0.015, 0.08, 'RMSE='+str(rmse_2000.crujra_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.015, 0.02, 'MAE='+str(mae_2000.crujra_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.115, 0.08, 'RMSE='+str(rmse_2000.crujra_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.115, 0.02, 'MAE='+str(mae_2000.crujra_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.22, 0.08, 'RMSE='+str(rmse_2000.crujra_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.22, 0.02, 'MAE='+str(mae_2000.crujra_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.32, 0.08, 'RMSE='+str(rmse_2000.crujraP_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.32, 0.02, 'MAE='+str(mae_2000.crujraP_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.42, 0.08, 'RMSE='+str(rmse_2000.crujraP_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.42, 0.02, 'MAE='+str(mae_2000.crujraP_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.52, 0.08, 'RMSE='+str(rmse_2000.crujraP_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.52, 0.02, 'MAE='+str(mae_2000.crujraP_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.625, 0.08, 'RMSE='+str(rmse_2000.oshd_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.625, 0.02, 'MAE='+str(mae_2000.oshd_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.725, 0.08, 'RMSE='+str(rmse_2000.oshd_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.725, 0.02, 'MAE='+str(mae_2000.oshd_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.82, 0.08, 'RMSE='+str(rmse_2000.oshd_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.82, 0.02, 'MAE='+str(mae_2000.oshd_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

plt.text(0.915, 0.08, 'RMSE='+str(rmse_2000.HS_jim.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.915, 0.02, 'MAE='+str(mae_2000.HS_jim.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)

axes = fig.add_subplot(313)
ax = sns.boxplot(data=(df_3000[['diff_cru_orig', 'diff_cru_new', 'diff_cru_nofor', 'diff_cruP_orig', 'diff_cruP_new',
                                'diff_cruP_nofor', 'diff_oshd_orig', 'diff_oshd_new', 'diff_oshd_nofor', 'FSM']]),
                 ax=axes, palette=my_pal, flierprops=flierprops, linewidth=1.6, saturation=0.8, showfliers=False,
                 meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": "6"},
                 medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), showmeans=True)
ax.axhline(0, color='k', linewidth=2.5)

axes.set_title(r"$\bf{(c)}$"+' Locations > 2000m (n = '+str(len(above_2000))+')', loc='left', fontsize=16)
axes.yaxis.grid(True)
axes.set_ylabel(u'Δ HS [m] ', fontsize=14)
for label in (axes.get_xticklabels() + axes.get_yticklabels()):
    label.set_fontsize(13)
axes.set_xticklabels([r'Clim$_{CRU^{pt}}$+LU$_{Gl}$', 'Clim$_{CRU^{pt}}$+LU$_{HR}$', 'Clim$_{CRU^{pt}}$+LU$_{nofor}$',
                      "Clim$_{CRU*^{pt}}$+LU$_{Gl}$", "Clim$_{CRU*^{pt}}$+LU$_{HR}$", "Clim$_{CRU*^{pt}}$+LU$_{nofor}$",
                      "Clim$_{OSHD^{pt}}$+LU$_{Gl}$", "Clim$_{OSHD^{pt}}$+LU$_{HR}$", "Clim$_{OSHD^{pt}}$+LU$_{nofor}$",
                      "Ref. (FSM2)"], rotation=20, ha='right')
axes.set_ylim([-2.1, 2.1])

for patch in ax.patches[1:-1:3]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.6))

for patch in ax.patches[2:-1:3]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0.2))

plt.text(0.015, 0.08, 'RMSE='+str(rmse_3000.crujra_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.015, 0.02, 'MAE='+str(mae_3000.crujra_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.115, 0.08, 'RMSE='+str(rmse_3000.crujra_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.115, 0.02, 'MAE='+str(mae_3000.crujra_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.22, 0.08, 'RMSE='+str(rmse_3000.crujra_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.22, 0.02, 'MAE='+str(mae_3000.crujra_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.32, 0.08, 'RMSE='+str(rmse_3000.crujraP_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.32, 0.02, 'MAE='+str(mae_3000.crujraP_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.42, 0.08, 'RMSE='+str(rmse_3000.crujraP_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.42, 0.02, 'MAE='+str(mae_3000.crujraP_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.52, 0.08, 'RMSE='+str(rmse_3000.crujraP_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.52, 0.02, 'MAE='+str(mae_3000.crujraP_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.625, 0.08, 'RMSE='+str(rmse_3000.oshd_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.625, 0.02, 'MAE='+str(mae_3000.oshd_origSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.725, 0.08, 'RMSE='+str(rmse_3000.oshd_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.725, 0.02, 'MAE='+str(mae_3000.oshd_newSurf.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.82, 0.08, 'RMSE='+str(rmse_3000.oshd_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.82, 0.02, 'MAE='+str(mae_3000.oshd_nofor.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.915, 0.08, 'RMSE='+str(rmse_3000.HS_jim.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.text(0.915, 0.02, 'MAE='+str(mae_3000.HS_jim.values[0])+'m', rotation=0, ha='left', va='bottom',
         transform=axes.transAxes, fontsize=8.3)
plt.tight_layout()
bf_out = bf / 'figures_final_revision'
fig.savefig(bf_out / 'f02.png', facecolor='white', transparent=False, bbox_inches='tight', dpi=1000)
fig.savefig(bf_out / 'f02.pdf', facecolor='white', transparent=False, bbox_inches='tight')

