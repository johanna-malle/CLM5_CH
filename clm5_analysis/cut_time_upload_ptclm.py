import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import dask
import netCDF4
import pandas as pd
import platform
import glob
import os


if platform.system() == 'Windows':
    bf_locs = Path('L:\malle\CLM5_CH')
else:
    bf_locs = Path('/home/lud11/malle/CLM5_CH')

path_FSM = bf_locs / 'FSM_new' / 'analysed_points'
all_files = glob.glob(os.path.join(path_FSM, "*.csv"))  # just do this once to get all ids

if platform.system() == 'Windows':
    all_locs_comp = list((f.split('\\')[-1]).split('_')[1] for f in all_files)
else:
    all_locs_comp = list((f.split('/')[-1]).split('_')[1] for f in all_files)

K = "MAE2"
K1 = "5DO"
all_locs = [i for i in all_locs_comp if (i != K and i != K1)]

bf = Path('C:/Users/malle/OneDrive - Eidg. Forschungsanstalt WSL/Documents/paper_clm5_ch/final_submission/'
          'upload_envidat')

bf_ptclm5 = bf / 'point_simulations_CLM5'

ptclm5_all = ['PTCLM_all_OSHD_newSurf', 'PTCLM_all_OSHD_origSurf', 'PTCLM_all_CRUJRA_Nolapse_newSurf',
              'PTCLM_all_CRUJRA_Nolapse_origSurf', 'PTCLM_all_CRUJRA_lapse_newSurf', 'PTCLM_all_CRUJRA_lapse_origSurf',
              'PTCLM5_nofor_CRU_lapse', 'PTCLM5_nofor_CRU_Nolapse', 'PTCLM5_nofor_OSHD']

name_out = ['PTCLM5_ClimOSHD_LUHR', 'PTCLM5_ClimOSHD_LUGl', 'PTCLM5_ClimCRU_LUHR',
            'PTCLM5_ClimCRU_LUGl', 'PTCLM5_ClimCRUstar_LUHR', 'PTCLM5_ClimCRUstar_LUGl',
            'PTCLM5_ClimCRUstar_LUnofor', 'PTCLM5_ClimCRU_LUnofor', 'PTCLM5_ClimOSHD_LUnofor']

for id_in, name_run in enumerate(ptclm5_all):
    print(name_run)
    bf_in = bf_ptclm5 / name_run
    bf_out = bf_ptclm5 / 'to_upload' / name_out[id_in]
    bf_out.mkdir(exist_ok=True, parents=True)

    for locs in all_locs:
        file_in = list(bf_in.glob(locs+'*.csv'))
        if file_in == []:
            file_in = list(bf_in.glob('*' + locs+'*.csv'))
        data_in = pd.read_csv(file_in[0], index_col='timeyears_clm')
        data_in.index = pd.to_datetime(data_in.index)
        filtered_df = data_in.loc[(data_in.index >= '2015-01-01') & (data_in.index <= '2020-01-01')]
        snow_related = filtered_df[['HS', 'SWE']]
        snow_related.rename(columns={'HS': 'SNOW_DEPTH [m]', 'SWE': 'SNOW_WATER_EQUIVALENT [mm]'}, inplace=True)
        snow_related.to_csv(bf_out / Path('PTCLM5_' + locs + '.csv'))
