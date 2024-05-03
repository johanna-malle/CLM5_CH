# -*- coding: utf-8 -*-
"""
Desc:
Created on 21.07.23 17:32
@author: malle
"""


import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import rioxarray
import seaborn as sns
import datetime
import glob
import pandas as pd
import platform


def open_file_all(run_in, id_veg_99, id_dem_3000, id_domain):
    print(run_in)
    if run_in.__contains__('deg') | run_in.__contains__('new'):
        gpp_in = xr.open_dataset(bf / run_in / 're_FPSN_SP_MONTH_SUM.nc')  # already in units of gC m-2 mo-1
        evap_in = xr.open_dataset(bf / run_in / 're_QFLX_EVAP_TOT.nc')
    else:
        gpp_in = xr.open_dataset(bf / run_in / 'FPSN_SP_MONTH_SUM.nc')  # already in units of gC m-2 mo-1
        evap_in = xr.open_dataset(glob.glob(str(bf / run_in) + '/QFLX_EVAP_TOT*')[0]) * 3600 * 24  # mm/s to mm/day

    evap_in_sum_mon = evap_in.resample(time='1M').sum()

    gpp_in.coords['mask_veg'] = (('lat', 'lon'), (id_veg_99 * 1).squeeze().data)
    gpp_in.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_dem_3000 * 1).squeeze().data))
    gpp_in.coords['domain'] = (('lat', 'lon'), (id_domain * 1).squeeze().data)
    datetimeindex = gpp_in.indexes['time'].to_datetimeindex()
    gpp_in['time'] = datetimeindex

    evap_in_sum_mon.coords['mask_veg'] = (('lat', 'lon'), (id_veg_99 * 1).squeeze().data)
    evap_in_sum_mon.coords['mask_3000'] = (('lat', 'lon'), np.flipud((id_dem_3000 * 1).squeeze().data))
    evap_in_sum_mon.coords['domain'] = (('lat', 'lon'), (id_domain * 1).squeeze().data)
    datetimeindex = evap_in_sum_mon.indexes['time'].to_datetimeindex()
    evap_in_sum_mon['time'] = datetimeindex

    gpp_avg = gpp_in.where((gpp_in.mask_veg == 1) & (gpp_in.domain == 0)).mean(dim=['lat', 'lon'], skipna=True)
    evap_avg = evap_in_sum_mon.where((evap_in_sum_mon.mask_veg == 1) &
                                     (evap_in_sum_mon.domain == 0)).mean(dim=['lat', 'lon'], skipna=True)

    gpp_avg_3000 = gpp_in.where((gpp_in.mask_veg == 1) & (gpp_in.domain == 0) &
                                (gpp_in.mask_3000 == 0)).mean(dim=['lat', 'lon'], skipna=True)
    evap_avg_3000 = evap_in_sum_mon.where((evap_in_sum_mon.mask_veg == 1) & (evap_in_sum_mon.mask_3000 == 0) &
                                          (evap_in_sum_mon.domain == 0)).mean(dim=['lat', 'lon'], skipna=True)

    return gpp_avg, evap_avg, gpp_avg_3000, evap_avg_3000


def make_dfs(oshd1km_hr, oshd1km_gl, oshd025_gl, oshd05_gl, cru1km_hr, cru1km_gl, cru025_gl, cru05_gl,
             cru_l_1km_hr, cru_l_1km_gl, cru_l_025_gl, cru_l_05_gl):
    df_oshd1km_hr = pd.Series(oshd1km_hr.DATA, name=r'Clim$_{OSHD 1km}$+LU$_{HR 1km}$', index=oshd1km_hr.time)
    df_oshd1km_gl = pd.Series(oshd1km_gl.DATA, name=r'Clim$_{OSHD 1km}$+LU$_{Gl 1km}$', index=oshd1km_gl.time)
    df_oshd025 = pd.Series(oshd025_gl.DATA, name=r'Clim$_{OSHD 0.25\degree}$+LU$_{Gl 0.25\degree}$',
                           index=oshd025_gl.time)
    df_oshd05 = pd.Series(oshd05_gl.DATA, name=r'Clim$_{OSHD 0.5\degree}$+LU$_{Gl 0.5\degree}$', index=oshd05_gl.time)
    df_cru1km_hr_l = pd.Series(cru_l_1km_hr.DATA, name=r'Clim$_{CRU* 1km}$+LU$_{HR 1km}$', index=cru_l_1km_hr.time)
    df_cru1km_gl_l = pd.Series(cru_l_1km_gl.DATA, name=r'Clim$_{CRU* 1km}$+LU$_{Gl 1km}$', index=cru_l_1km_gl.time)
    df_cru025_l = pd.Series(cru_l_025_gl.DATA, name=r'Clim$_{CRU* 0.25\degree}$+LU$_{Gl 0.25\degree}$',
                            index=cru_l_025_gl.time)
    df_crud05_l = pd.Series(cru_l_05_gl.DATA, name=r'Clim$_{CRU* 0.5\degree}$+LU$_{Gl 0.5\degree}$',
                            index=cru_l_05_gl.time)
    df_cru1km_hr = pd.Series(cru1km_hr.DATA, name=r'Clim$_{CRU 1km}$+LU$_{HR 1km}$', index=cru1km_hr.time)
    df_cru1km_gl = pd.Series(cru1km_gl.DATA, name=r'Clim$_{CRU 1km}$+LU$_{Gl 1km}$', index=cru1km_gl.time)
    df_cru025 = pd.Series(cru025_gl.DATA, name=r'Clim$_{CRU 0.25\degree}$+LU$_{Gl 0.25\degree}$', index=cru025_gl.time)
    df_cru05 = pd.Series(cru05_gl.DATA, name=r'Clim$_{CRU 0.5\degree}$+LU$_{Gl 0.5\degree}$', index=cru05_gl.time)
    df_all = pd.concat([df_oshd1km_hr, df_oshd1km_gl, df_oshd025, df_oshd05,
                        df_cru1km_hr_l, df_cru1km_gl_l, df_cru025_l, df_crud05_l,
                        df_cru1km_hr, df_cru1km_gl, df_cru025, df_cru05], axis=1)
    return df_all


if __name__ == '__main__':

    if platform.system() == 'Windows':
        bf = Path('L:\malle\CLM5_CH')
    else:
        bf = Path('/home/lud11/malle/CLM5_CH')

    bf_plots = bf / 'snow_eco_new/gpp_snow/'

    surf_in = xr.open_dataset(bf / 'surfdata_1km_CH_v3_hist_16pfts_Irrig_CMIP6_NEW.nc')
    id_surf_99 = surf_in.PCT_NATVEG > 99

    oshd_in = xr.open_dataset(bf / 'OSHD_FILES' / 'FPSN_SP_MONTH_SUM.nc')
    id_dom = np.isnan(oshd_in.DATA.isel(time=0))

    dem = rioxarray.open_rasterio(bf / 'BAFU_DEM_2020_1000.tif')
    id_1000 = (dem > 0) & (dem < 1000)
    id_2000 = (dem >= 1000) & (dem < 2000)
    id_3000 = (dem >= 2000)

    my_pal = {r'Clim$_{CRU 1km}$+LU$_{HR 1km}$': (217 / 255, 95 / 255, 2 / 255),
              r'Clim$_{CRU 1km}$+LU$_{Gl 1km}$': (217 / 255, 95 / 255, 2 / 255),
              r'Clim$_{CRU 0.25\degree}$+LU$_{Gl 0.25\degree}$': (217 / 255, 95 / 255, 2 / 255),
              r'Clim$_{CRU 0.5\degree}$+LU$_{Gl 0.5\degree}$': (217 / 255, 95 / 255, 2 / 255),

              r'Clim$_{CRU* 1km}$+LU$_{HR 1km}$': (117 / 255, 112 / 255, 179 / 255),
              r'Clim$_{CRU* 1km}$+LU$_{Gl 1km}$': (117 / 255, 112 / 255, 179 / 255),
              r'Clim$_{CRU* 0.25\degree}$+LU$_{Gl 0.25\degree}$': (117 / 255, 112 / 255, 179 / 255),
              r'Clim$_{CRU* 0.5\degree}$+LU$_{Gl 0.5\degree}$': (117 / 255, 112 / 255, 179 / 255),

              r'Clim$_{OSHD 1km}$+LU$_{HR 1km}$': (27 / 255, 158 / 255, 119 / 255),
              r'Clim$_{OSHD 1km}$+LU$_{Gl 1km}$': (27 / 255, 158 / 255, 119 / 255),
              r'Clim$_{OSHD 0.25\degree}$+LU$_{Gl 0.25\degree}$': (27 / 255, 158 / 255, 119 / 255),
              r'Clim$_{OSHD 0.5\degree}$+LU$_{Gl 0.5\degree}$': (27 / 255, 158 / 255, 119 / 255)}

    gpp_oshd1km_HR, evap_oshd1km_HR, gpp_oshd1km_HR_3000, evap_oshd1km_HR_3000 = \
        open_file_all('OSHD_FILES', id_surf_99, id_3000, id_dom)
    gpp_oshd1km_GL, evap_oshd1km_GL, gpp_oshd1km_GL_3000, evap_oshd1km_GL_3000 = \
        open_file_all('OSHD_FILES_OLD', id_surf_99, id_3000, id_dom)
    gpp_oshd025, evap_oshd025, gpp_oshd025_3000, evap_oshd025_3000 = \
        open_file_all('OSHD_FILES_025_new', id_surf_99, id_3000, id_dom)
    gpp_oshd05, evap_oshd05, gpp_oshd05_3000, evap_oshd05_3000 = \
        open_file_all('OSHD_FILES_05_new', id_surf_99, id_3000, id_dom)

    gpp_cruL1km_HR, evap_cruL1km_HR, gpp_cruL1km_HR_3000, evap_cruL1km_HR_3000 = \
        open_file_all('CRUJRA_FILES', id_surf_99, id_3000, id_dom)
    gpp_cruL1km_GL, evap_cruL1km_GL, gpp_cruL1km_GL_3000, evap_cruL1km_GL_3000 = \
        open_file_all('CRUJRA_FILES_OLD', id_surf_99, id_3000, id_dom)
    gpp_cruL025, evap_cruL025, gpp_cruL025_3000, evap_cruL025_3000 = \
        open_file_all('CRUJRA_FILES_025deg_cru_new_lapse', id_surf_99, id_3000, id_dom)
    gpp_cruL05, evap_cruL05, gpp_cruL05_3000, evap_cruL05_3000 = \
        open_file_all('CRUJRA_FILES_05deg_cru_new_lapse', id_surf_99, id_3000, id_dom)

    gpp_cru1km_HR, evap_cru1km_HR, gpp_cru1km_HR_3000, evap_cru1km_HR_3000 = \
        open_file_all('CRUJRA_FILES_noLapse', id_surf_99, id_3000, id_dom)
    gpp_cru1km_GL, evap_cru1km_GL, gpp_cru1km_GL_3000, evap_cru1km_GL_3000 = \
        open_file_all('CRUJRA_FILES_noLapse_OLD', id_surf_99, id_3000, id_dom)
    gpp_cru025, evap_cru025, gpp_cru025_3000, evap_cru025_3000 = \
        open_file_all('CRUJRA_FILES_025deg_cru_new', id_surf_99, id_3000, id_dom)
    gpp_cru05, evap_cru05, gpp_cru05_3000, evap_cru05_3000 = \
        open_file_all('CRUJRA_FILES_05deg_cru_new', id_surf_99, id_3000, id_dom)

    df_gpp = make_dfs(gpp_oshd1km_HR, gpp_oshd1km_GL, gpp_oshd025, gpp_oshd05, gpp_cru1km_HR, gpp_cru1km_GL, gpp_cru025,
                      gpp_cru05, gpp_cruL1km_HR, gpp_cruL1km_GL, gpp_cruL025, gpp_cruL05)

    df_evap = make_dfs(evap_oshd1km_HR, evap_oshd1km_GL, evap_oshd025, evap_oshd05, evap_cru1km_HR, evap_cru1km_GL,
                       evap_cru025, evap_cru05, evap_cruL1km_HR, evap_cruL1km_GL, evap_cruL025, evap_cruL05)

    df_gpp_3000 = make_dfs(gpp_oshd1km_HR_3000, gpp_oshd1km_GL_3000, gpp_oshd025_3000, gpp_oshd05_3000,
                           gpp_cru1km_HR_3000, gpp_cru1km_GL_3000, gpp_cru025_3000, gpp_cru05_3000, gpp_cruL1km_HR_3000,
                           gpp_cruL1km_GL_3000, gpp_cruL025_3000, gpp_cruL05_3000)

    df_evap_3000 = make_dfs(evap_oshd1km_HR_3000, evap_oshd1km_GL_3000, evap_oshd025_3000, evap_oshd05_3000,
                            evap_cru1km_HR_3000, evap_cru1km_GL_3000, evap_cru025_3000, evap_cru05_3000,
                            evap_cruL1km_HR_3000, evap_cruL1km_GL_3000, evap_cruL025_3000, evap_cruL05_3000)

    dstart = datetime.datetime(2017, 1, 1)
    dend = datetime.datetime(2020, 1, 1)
    dstart_zoom = datetime.datetime(2018, 5, 15)
    dend_zoom = datetime.datetime(2018, 8, 15)

    fig = plt.figure(figsize=(10, 7))
    axes = fig.add_subplot(211)
    spl = sns.lineplot(data=df_gpp, ax=axes, palette=my_pal)
    axes.set_xlim([dstart, dend])
    axes.yaxis.grid(True)
    axes.xaxis.grid(True)
    axes.set_ylabel('GPP [gC m$^{-2}$ month$^{-1}$]')
    for line in spl.lines[::4]:
        line.set_linestyle("-")
    for line in spl.lines[1::4]:
        line.set_linestyle("--")
    for line in spl.lines[2::4]:
        line.set_linestyle("-.")
    for line in spl.lines[3::4]:
        line.set_linestyle(":")
    spl.axvline(dstart_zoom, linestyle='--', color='cornflowerblue', linewidth=0.9)
    spl.axvline(dend_zoom, linestyle='--', color='cornflowerblue', linewidth=0.9)
    spl.legend(loc='upper left', bbox_to_anchor=(0.1, 1.45), ncol=3)

    axes = fig.add_subplot(212)
    spl = sns.lineplot(data=df_evap, ax=axes, palette=my_pal, legend=False)
    axes.set_xlim([dstart, dend])
    axes.yaxis.grid(True)
    axes.xaxis.grid(True)
    axes.set_ylabel('ET [mm month$^{-1}$]')
    for line in spl.lines[::4]:
        line.set_linestyle("-")
    for line in spl.lines[1::4]:
        line.set_linestyle("--")
    for line in spl.lines[2::4]:
        line.set_linestyle("-.")
    for line in spl.lines[3::4]:
        line.set_linestyle(":")
    plt.tight_layout()
    spl.axvline(dstart_zoom, linestyle='--', color='cornflowerblue', linewidth=0.9)
    spl.axvline(dend_zoom, linestyle='--', color='cornflowerblue', linewidth=0.9)
    plt.show()
    fig.savefig(bf_plots / 'temporal_gpp_evap.png', dpi=500, bbox_inches='tight', facecolor='white', transparent=False)

    fig = plt.figure(figsize=(10, 7))
    axes = fig.add_subplot(211)
    spl = sns.lineplot(data=df_gpp_3000, ax=axes, palette=my_pal)
    axes.set_xlim([dstart, dend])
    axes.yaxis.grid(True)
    axes.xaxis.grid(True)
    for line in spl.lines[::4]:
        line.set_linestyle("-")
    for line in spl.lines[1::4]:
        line.set_linestyle("--")
    for line in spl.lines[2::4]:
        line.set_linestyle("-.")
    for line in spl.lines[3::4]:
        line.set_linestyle(":")
    spl.legend(loc='upper left', bbox_to_anchor=(0.1, 1.45), ncol=3)
    spl.axvline(dstart_zoom, linestyle='--', color='cornflowerblue', linewidth=0.9)
    spl.axvline(dend_zoom, linestyle='--', color='cornflowerblue', linewidth=0.9)
    axes.set_ylabel('GPP [gC m$^{-2}$ month$^{-1}$]')

    axes = fig.add_subplot(212)
    spl = sns.lineplot(data=df_evap_3000, ax=axes, palette=my_pal, legend=False)
    axes.set_xlim([dstart, dend])
    axes.yaxis.grid(True)
    axes.xaxis.grid(True)
    axes.set_ylabel('ET [mm month$^{-1}$]')
    for line in spl.lines[::4]:
        line.set_linestyle("-")
    for line in spl.lines[1::4]:
        line.set_linestyle("--")
    for line in spl.lines[2::4]:
        line.set_linestyle("-.")
    for line in spl.lines[3::4]:
        line.set_linestyle(":")
    plt.tight_layout()
    spl.axvline(dstart_zoom, linestyle='--', color='cornflowerblue', linewidth=0.9)
    spl.axvline(dend_zoom, linestyle='--', color='cornflowerblue', linewidth=0.9)
    plt.show()
    fig.savefig(bf_plots / 'temporal_gpp_evap_3000.png', dpi=500, bbox_inches='tight', facecolor='white',
                transparent=False)

    # zoom plot
    fig = plt.figure(figsize=(5, 7))
    axes = fig.add_subplot(211)
    spl = sns.lineplot(data=df_gpp, ax=axes, palette=my_pal, legend=False)
    axes.set_xlim([dstart_zoom, dend_zoom])
    axes.set_ylim([150, 240])
    axes.yaxis.grid(True)
    axes.xaxis.grid(True)
    axes.set_ylabel('GPP [gC m$^{-2}$ month$^{-1}$]')
    for line in spl.lines[::4]:
        line.set_linestyle("-")
    for line in spl.lines[1::4]:
        line.set_linestyle("--")
    for line in spl.lines[2::4]:
        line.set_linestyle("-.")
    for line in spl.lines[3::4]:
        line.set_linestyle(":")
    axes.tick_params(labelbottom=False)

    axes = fig.add_subplot(212)
    spl = sns.lineplot(data=df_evap, ax=axes, palette=my_pal, legend=False)
    axes.set_xlim([dstart_zoom, dend_zoom])
    axes.set_ylim([30, 130])
    axes.yaxis.grid(True)
    axes.xaxis.grid(True)
    axes.set_ylabel('ET [mm month$^{-1}$]')
    for line in spl.lines[::4]:
        line.set_linestyle("-")
    for line in spl.lines[1::4]:
        line.set_linestyle("--")
    for line in spl.lines[2::4]:
        line.set_linestyle("-.")
    for line in spl.lines[3::4]:
        line.set_linestyle(":")
    axes.tick_params(axis='x', rotation=90)
    plt.tight_layout()
    plt.show()
    fig.savefig(bf_plots / 'temporal_gpp_evap_zoom.png', dpi=500, bbox_inches='tight', facecolor='white',
                transparent=False)

    fig = plt.figure(figsize=(5, 7))
    axes = fig.add_subplot(211)
    spl = sns.lineplot(data=df_gpp_3000, ax=axes, palette=my_pal, legend=False)
    axes.set_xlim([dstart_zoom, dend_zoom])
    axes.set_ylim([150, 250])
    axes.yaxis.grid(True)
    axes.xaxis.grid(True)
    axes.set_ylabel('GPP [gC m$^{-2}$ month$^{-1}$]')
    for line in spl.lines[::4]:
        line.set_linestyle("-")
    for line in spl.lines[1::4]:
        line.set_linestyle("--")
    for line in spl.lines[2::4]:
        line.set_linestyle("-.")
    for line in spl.lines[3::4]:
        line.set_linestyle(":")
    axes.tick_params(labelbottom=False)

    axes = fig.add_subplot(212)
    spl = sns.lineplot(data=df_evap_3000, ax=axes, palette=my_pal, legend=False)
    axes.set_xlim([dstart_zoom, dend_zoom])
    axes.set_ylim([30, 130])
    axes.yaxis.grid(True)
    axes.xaxis.grid(True)
    axes.set_ylabel('ET [mm month$^{-1}$]')
    for line in spl.lines[::4]:
        line.set_linestyle("-")
    for line in spl.lines[1::4]:
        line.set_linestyle("--")
    for line in spl.lines[2::4]:
        line.set_linestyle("-.")
    for line in spl.lines[3::4]:
        line.set_linestyle(":")
    axes.tick_params(axis='x', rotation=90)
    plt.tight_layout()
    plt.show()
    fig.savefig(bf_plots / 'temporal_gpp_evap_zoom_3000.png', dpi=500, bbox_inches='tight', facecolor='white',
                transparent=False)
