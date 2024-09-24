import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import dask
import netCDF4


bf = Path('C:/Users/malle/OneDrive - Eidg. Forschungsanstalt WSL/Documents/paper_clm5_ch/final_submission/'
          'upload_envidat')

bf_gridded = bf / 'gridded_simulations_CLM5'
bf_ptclm5 = bf / 'point_simulations_CLM5'

run_in_all = ['CRUJRA_FILES_noLapse_OLD', 'CRUJRA_FILES_noLapse', 'CRUJRA_FILES_OLD', 'CRUJRA_FILES',
              'OSHD_FILES_OLD', 'OSHD_FILES',
              'CRUJRA_FILES_05deg_cru_new', 'CRUJRA_FILES_05deg_cru_new_lapse', 'OSHD_FILES_05_new',
              'CRUJRA_FILES_025deg_cru_new', 'CRUJRA_FILES_025deg_cru_new_lapse', 'OSHD_FILES_025_new']

name_out = ['ClimCRU_LUGl_1km', 'ClimCRU_LUHR_1km', 'ClimCRUstar_LUGl_1km', 'ClimCRUstar_LUHR_1km',
            'ClimOSHD_LUGl_1km', 'ClimOSHD_LUHR_1km',
            'ClimCRU_LUGl_05deg', 'ClimCRUstar_LUGl_05deg', 'ClimOSHD_LUGl_05deg',
            'ClimCRU_LUGl_025deg', 'ClimCRUstar_LUGl_025deg', 'ClimOSHD_LUGl_025deg']

name_in_all_snow = ['SNOW_DEPTH_SP_CRUJRA_noLapse.nc', 'SNOW_DEPTH_SP_CRUJRA_noLapse.nc',
                    'SNOW_DEPTH_SP_CRUJRA.nc', 'SNOW_DEPTH_SP_CRUJRA.nc',
                    'SNOW_DEPTH_SP_OSHD.nc', 'SNOW_DEPTH_SP_OSHD.nc',
                    'SNOW_DEPTH.nc', 'SNOW_DEPTH.nc', 'SNOW_DEPTH.nc', 'SNOW_DEPTH.nc',
                    'SNOW_DEPTH.nc', 'SNOW_DEPTH.nc']

name_in_all_qflx = ['QFLX_EVAP_TOT_SP_MONTH_SUM.nc', 'QFLX_EVAP_TOT_SP_MONTH_SUM.nc',
                    'QFLX_EVAP_TOT_SP_MONTH_SUM.nc', 'QFLX_EVAP_TOT_SP_MONTH_SUM.nc',
                    'QFLX_EVAP_TOT_SP_MONTH_SUM.nc', 'QFLX_EVAP_TOT_SP_MONTH_SUM.nc',
                    'QFLX_EVAP_TOT.nc', 'QFLX_EVAP_TOT.nc', 'QFLX_EVAP_TOT.nc', 'QFLX_EVAP_TOT.nc', 'QFLX_EVAP_TOT.nc',
                    'QFLX_EVAP_TOT.nc']

name_in_all_fpsn = ['FPSN_SP_MONTH_SUM.nc', 'FPSN_SP_MONTH_SUM.nc',
                    'FPSN_SP_MONTH_SUM.nc', 'FPSN_SP_MONTH_SUM.nc',
                    'FPSN_SP_MONTH_SUM.nc', 'FPSN_SP_MONTH_SUM.nc',
                    'FPSN.nc', 'FPSN.nc', 'FPSN.nc', 'FPSN.nc', 'FPSN.nc', 'FPSN.nc']

for id_in, name_run in enumerate(run_in_all):
    print(name_run)
    bf_in = bf_gridded / name_run
    bf_out = bf_gridded / 'to_upload' / name_out[id_in]
    bf_out.mkdir(exist_ok=True, parents=True)
    name_in_snow = name_in_all_snow[id_in]
    name_in_fpsn = name_in_all_fpsn[id_in]
    name_in_qflx = name_in_all_qflx[id_in]

    snow_in = xr.open_dataset(bf_in / name_in_snow)
    snow_cut = snow_in.isel(time=slice(1825, np.shape(snow_in.time)[0]))
    snow_cut.to_netcdf(bf_out / 'SNOW_DEPTH.nc')

    fpsn = xr.open_dataset(bf_in / name_in_fpsn)
    fpsn_cut = fpsn.isel(time=slice(60, np.shape(fpsn.time)[0]))
    fpsn_cut.to_netcdf(bf_out / 'FPSN_MONTH_SUM.nc')

    qflx = xr.open_dataset(bf_in / name_in_qflx)
    qflx_cut = qflx.isel(time=slice(60, np.shape(qflx.time)[0]))
    qflx_cut.to_netcdf(bf_out / 'QFLX_EVAP_TOT_MONTH_SUM.nc')