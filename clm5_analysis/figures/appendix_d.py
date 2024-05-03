# -*- coding: utf-8 -*-
"""
Desc:
Created on 13.07.23 15:34
@author: malle
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib
import platform
import seaborn as sns
import pandas as pd


def make_rel_plots(oshd1km_hr, oshd1km_gl, oshd025_gl, oshd05_gl, cru1km_hr, cru1km_gl, cru025_gl, cru05_gl,
                   crulap1km_hr, crulap1km_gl, crulap025_gl, crulap05_gl, units_in, save_in, plot_dir):
    proj_data = ccrs.PlateCarree()
    proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)

    vmin = np.nanmin([(oshd1km_gl.DATA - oshd1km_hr.DATA),
                      (oshd025_gl.DATA.data - oshd1km_hr.DATA), (oshd05_gl.DATA.data - oshd1km_hr.DATA),
                      (cru1km_gl.DATA - oshd1km_hr.DATA), (cru1km_hr.DATA - oshd1km_hr.DATA),
                      (crulap1km_gl.DATA - oshd1km_hr.DATA), (crulap1km_hr.DATA - oshd1km_hr.DATA),
                      (cru025_gl.DATA.data - oshd1km_hr.DATA), (cru05_gl.DATA.data - oshd1km_hr.DATA),
                      (crulap05_gl.DATA.data - oshd1km_hr.DATA), (crulap025_gl.DATA.data - oshd1km_hr.DATA)])

    vmax = np.nanmax([(oshd1km_gl.DATA - oshd1km_hr.DATA),
                      (oshd025_gl.DATA.data - oshd1km_hr.DATA), (oshd05_gl.DATA.data - oshd1km_hr.DATA),
                      (cru1km_gl.DATA - oshd1km_hr.DATA), (cru1km_hr.DATA - oshd1km_hr.DATA),
                      (crulap1km_gl.DATA - oshd1km_hr.DATA), (crulap1km_hr.DATA - oshd1km_hr.DATA),
                      (cru025_gl.DATA.data - oshd1km_hr.DATA), (cru05_gl.DATA.data - oshd1km_hr.DATA),
                      (crulap05_gl.DATA.data - oshd1km_hr.DATA), (crulap025_gl.DATA.data - oshd1km_hr.DATA)])
    vabs = np.max(np.abs((vmin, vmax)))

    diff_copy = oshd1km_hr.copy()

    fig, axs = plt.subplots(3, 4, figsize=[13, 11], layout='constrained', subplot_kw={'projection': proj_map})
    plt.suptitle('Year: ' + str(yr_in), fontweight='bold', fontsize=15)

    a1 = oshd1km_hr.DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[0, 0])
    axs[0, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[0, 0].axis('off')
    axs[0, 0].set_title(r'Ref.: Clim$\bf{_{OSHD 1km}}$+LU$\bf{_{HR 1km}}$', fontweight='bold', fontsize=14)
    cbar = fig.colorbar(a1, ax=axs[0, 0], shrink=0.6, location='bottom')
    cbar.set_label(units_in, rotation=0, labelpad=5)

    (oshd1km_gl.DATA - oshd1km_hr.DATA).plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False,
                                             vmin=-vabs, vmax=vabs, ax=axs[0, 1])
    axs[0, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[0, 1].axis('off')
    axs[0, 1].set_title('Clim$_{OSHD 1km}$+LU$_{Gl 1km}$ - Ref.')

    diff_copy.DATA.data = (oshd025_gl.DATA.data - oshd1km_hr.DATA.data)
    diff_copy.DATA.plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False, vmin=-vabs, vmax=vabs, ax=axs[0, 2])
    axs[0, 2].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[0, 2].axis('off')
    axs[0, 2].set_title(r'Clim$_{OSHD 0.25\degree}$+LU$_{Gl 0.25\degree}$ - Ref.')

    diff_copy.DATA.data = (oshd05_gl.DATA.data - oshd1km_hr.DATA.data)
    diff_copy.DATA.plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False, vmin=-vabs, vmax=vabs, ax=axs[0, 3])
    axs[0, 3].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[0, 3].axis('off')
    axs[0, 3].set_title(r'Clim$_{OSHD 0.5\degree}$+LU$_{Gl 0.5\degree}$ - Ref.')

    # now cru lapse
    (crulap1km_hr.DATA - oshd1km_hr.DATA).plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False,
                                               vmin=-vabs, vmax=vabs, ax=axs[1, 0])
    axs[1, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[1, 0].axis('off')
    axs[1, 0].set_title('Clim$_{CRU* 1km}$+LU$_{HR 1km}$ - Ref.')

    (crulap1km_gl.DATA - oshd1km_hr.DATA).plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False,
                                               vmin=-vabs, vmax=vabs, ax=axs[1, 1])
    axs[1, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[1, 1].axis('off')
    axs[1, 1].set_title('Clim$_{CRU* 1km}$+LU$_{Gl 1km}$ - Ref.')

    diff_copy.DATA.data = (crulap025_gl.DATA.data - oshd1km_hr.DATA.data)
    diff_copy.DATA.plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False, vmin=-vabs, vmax=vabs, ax=axs[1, 2])
    axs[1, 2].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[1, 2].axis('off')
    axs[1, 2].set_title(r'Clim$_{CRU* 0.25\degree}$+LU$_{Gl 0.25\degree}$ - Ref.')

    diff_copy.DATA.data = (crulap05_gl.DATA.data - oshd1km_hr.DATA.data)
    diff_copy.DATA.plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False, vmin=-vabs, vmax=vabs, ax=axs[1, 3])
    axs[1, 3].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[1, 3].axis('off')
    axs[1, 3].set_title(r'Clim$_{CRU* 0.5\degree}$+LU$_{Gl 0.5\degree}$ - Ref.')

    # now cru
    (cru1km_hr.DATA - oshd1km_hr.DATA).plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False,
                                            vmin=-vabs, vmax=vabs, ax=axs[2, 0])
    axs[2, 0].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[2, 0].axis('off')
    axs[2, 0].set_title('Clim$_{CRU 1km}$+LU$_{HR 1km}$ - Ref.')

    (cru1km_gl.DATA - oshd1km_hr.DATA).plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False,
                                            vmin=-vabs, vmax=vabs, ax=axs[2, 1])
    axs[2, 1].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[2, 1].axis('off')
    axs[2, 1].set_title('Clim$_{CRU 1km}$+LU$_{Gl 1km}$ - Ref.')

    diff_copy.DATA.data = (cru025_gl.DATA.data - oshd1km_hr.DATA.data)
    diff_copy.DATA.plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False, vmin=-vabs, vmax=vabs, ax=axs[2, 2])
    axs[2, 2].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[2, 2].axis('off')
    axs[2, 2].set_title(r'Clim$_{CRU 0.25\degree}$+LU$_{Gl 0.25\degree}$ - Ref.')

    diff_copy.DATA.data = (cru05_gl.DATA.data - oshd1km_hr.DATA.data)
    a2 = diff_copy.DATA.plot(transform=proj_data, cmap='RdBu_r', add_colorbar=False, vmin=-vabs, vmax=vabs,
                             ax=axs[2, 3])
    axs[2, 3].add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
    axs[2, 3].axis('off')
    axs[2, 3].set_title(r'Clim$_{CRU 0.5\degree}$+LU$_{Gl 0.5\degree}$ - Ref.')
    cbar = fig.colorbar(a2, ax=axs[2, 0:4], shrink=0.75, location='bottom', aspect=30)
    cbar.set_label('$\Delta$ ' + units_in, rotation=0, labelpad=5)

    fig.savefig(plot_dir / Path(save_in + str(yr_in) + '.png'),
                facecolor='white', transparent=False)
    plt.close()


def make_df_boxplots(oshd1km_hr, oshd1km_gl, oshd025_gl, oshd05_gl, cru1km_hr, cru1km_gl, cru025_gl, cru05_gl,
                     crulap1km_hr, crulap1km_gl, crulap025_gl, crulap05_gl):
    df_oshd1kmhr = pd.Series(oshd1km_hr.where((oshd1km_hr.mask_veg == 1)).DATA.data.flatten(),
                             name=r'Clim$_{OSHD 1km}$+LU$_{HR 1km}$')
    df_oshd1kmgl = pd.Series(oshd1km_gl.where((oshd1km_gl.mask_veg == 1)).DATA.data.flatten(),
                             name=r'Clim$_{OSHD 1km}$+LU$_{Gl 1km}$')
    df_oshd025 = pd.Series(oshd025_gl.where((oshd025_gl.mask_veg == 1)).DATA.data.flatten(),
                           name=r'Clim$_{OSHD 0.25\degree}$+LU$_{Gl 0.25\degree}$')
    df_oshd05 = pd.Series(oshd05_gl.where((oshd05_gl.mask_veg == 1)).DATA.data.flatten(),
                          name=r'Clim$_{OSHD 0.5\degree}$+LU$_{Gl 0.5\degree}$')
    df_cru1kmhr_lap = pd.Series(crulap1km_hr.where((crulap1km_hr.mask_veg == 1)).DATA.data.flatten(),
                                name=r'Clim$_{CRU* 1km}$+LU$_{HR 1km}$')
    df_cru1kmgl_lap = pd.Series(crulap1km_gl.where((crulap1km_gl.mask_veg == 1)).DATA.data.flatten(),
                                name=r'Clim$_{CRU* 1km}$+LU$_{Gl 1km}$')
    df_cru025_lap = pd.Series(crulap025_gl.where((crulap025_gl.mask_veg == 1)).DATA.data.flatten(),
                              name=r'Clim$_{CRU* 0.25\degree}$+LU$_{Gl 0.25\degree}$')
    df_crud05_lap = pd.Series(crulap05_gl.where((crulap05_gl.mask_veg == 1)).DATA.data.flatten(),
                              name=r'Clim$_{CRU* 0.5\degree}$+LU$_{Gl 0.5\degree}$')
    df_cru1kmhr = pd.Series(cru1km_hr.where((cru1km_hr.mask_veg == 1)).DATA.data.flatten(),
                            name=r'Clim$_{CRU 1km}$+LU$_{HR 1km}$')
    df_cru1kmgl = pd.Series(cru1km_gl.where((cru1km_gl.mask_veg == 1)).DATA.data.flatten(),
                            name=r'Clim$_{CRU 1km}$+LU$_{Gl 1km}$')
    df_cru025 = pd.Series(cru025_gl.where((cru025_gl.mask_veg == 1)).DATA.data.flatten(),
                          name=r'Clim$_{CRU 0.25\degree}$+LU$_{Gl 0.25\degree}$')
    df_cru05 = pd.Series(cru05_gl.where((cru05_gl.mask_veg == 1)).DATA.data.flatten(),
                         name=r'Clim$_{CRU 0.5\degree}$+LU$_{Gl 0.5\degree}$')
    df_all = pd.concat([df_oshd1kmhr, df_oshd1kmgl, df_oshd025, df_oshd05,
                        df_cru1kmhr_lap, df_cru1kmgl_lap, df_cru025_lap, df_crud05_lap,
                        df_cru1kmhr, df_cru1kmgl, df_cru025, df_cru05], axis=1)
    df_all.loc[(np.isnan(df_all[r'Clim$_{OSHD 1km}$+LU$_{HR 1km}$']))] = np.nan
    return df_all


def open_files(run_in, bf, year_in):
    ending = str(year_in) + '_v0.nc'
    file_in_gpp = 'gpp_avg_' + ending
    file_in_gpp_0506 = 'gpp_avg_0506_' + ending
    file_in_snow = 'snow_count_' + ending
    file_in_snow_yr = 'snow_count_year_' + ending
    file_in_snow_hydroyr = 'snow_count_HYDROyear_' + ending
    file_in_evap_yr = 'total_ET_year_' + ending
    file_in_evap_hydroyr = 'total_ET_HYDORyear_' + ending
    file_in_snow_total = 'snow_total_' + ending
    file_in_snow_total_yr = 'snow_total_year_' + ending
    file_in_melt_date = 'melt_out_' + ending
    file_in_melt_doy = 'melt_doy_' + ending
    file_in_swe_total = 'swe_inc_total_HYDROyear_' + ending

    gpp = xr.open_dataset(bf / run_in / file_in_gpp)
    gpp0506 = xr.open_dataset(bf / run_in / file_in_gpp_0506)
    snow = xr.open_dataset(bf / run_in / file_in_snow)
    snow_yr = xr.open_dataset(bf / run_in / file_in_snow_yr)
    snow_hyr = xr.open_dataset(bf / run_in / file_in_snow_hydroyr)
    evap = xr.open_dataset(bf / run_in / file_in_evap_yr)
    evap_hyr = xr.open_dataset(bf / run_in / file_in_evap_hydroyr)
    snow_total = xr.open_dataset(bf / run_in / file_in_snow_total)
    snow_total_yr = xr.open_dataset(bf / run_in / file_in_snow_total_yr)
    swe_total_hyr = xr.open_dataset(bf / run_in / file_in_swe_total)
    melt_out = xr.open_dataset(bf / run_in / file_in_melt_date)
    melt_doy = xr.open_dataset(bf / run_in / file_in_melt_doy)

    # set melt_out dates to nan when < 50 snow days per hydrological year or when later than July 1st
    if run_in == 'OSHD_FILES':
        melt_out = melt_out.where(snow_hyr.DATA > 50, np.datetime64('NaT'))
        melt_doy = melt_doy.where(snow_hyr.DATA > 50, np.nan)
        melt_doy = melt_doy.where(melt_doy.DATA < 183, np.nan)
        melt_doy = melt_doy.where(melt_doy.DATA > 31, np.nan)

    return gpp, gpp0506, evap, evap_hyr, melt_out, melt_doy, \
        snow, snow_yr, snow_hyr, snow_total, snow_total_yr, swe_total_hyr


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

if __name__ == '__main__':
    # make switch to windows if working from home
    if platform.system() == 'Windows':
        bf = Path('L:\malle\CLM5_CH')
    else:
        bf = Path('/home/lud11/malle/CLM5_CH')

    bf_plots = bf / 'snow_eco_new/gpp_snow/relative_spatial_plots'

    yr_all = [2016, 2017, 2018, 2019]

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

    for yr_in in yr_all:
        print(yr_in)
        gpp_oshd1km_HR, gpp0506_oshd1km_HR, evap_oshd1km_HR, evap_oshd1km_HR_HY, melt_out_oshd1km_HR, melt_doy_oshd1km_HR, \
            snow_oshd1km_HR, snow_yr_oshd1km_HR, snow_Hyr_oshd1km_HR, snow_total_oshd1km_HR, snow_total_yr_oshd1km_HR, \
            swe_total_Hyr_oshd1km_HR = open_files('OSHD_FILES', bf, yr_in)
        gpp_oshd1km_GL, gpp0506_oshd1km_GL, evap_oshd1km_GL, evap_oshd1km_GL_HY, melt_out_oshd1km_GL, melt_doy_oshd1km_GL, \
            snow_oshd1km_GL, snow_yr_oshd1km_GL, snow_Hyr_oshd1km_GL, snow_total_oshd1km_GL, snow_total_yr_oshd1km_GL, \
            swe_total_Hyr_oshd1km_GL = open_files('OSHD_FILES_OLD', bf, yr_in)
        gpp_oshd05_GL, gpp0506_oshd05_GL, evap_oshd05_GL, evap_oshd05_GL_HY, melt_out_oshd05_GL, melt_doy_oshd05_GL, \
            snow_oshd05_GL, snow_yr_oshd05_GL, snow_Hyr_oshd05_GL, snow_total_oshd05_GL, snow_total_yr_oshd05_GL, \
            swe_total_Hyr_oshd05_GL = open_files('OSHD_FILES_05_new', bf, yr_in)
        gpp_oshd025_GL, gpp0506_oshd025_GL, evap_oshd025_GL, evap_oshd025_GL_HY, melt_out_oshd025_GL, melt_doy_oshd025_GL, \
            snow_oshd025_GL, snow_yr_oshd025_GL, snow_Hyr_oshd025_GL, snow_total_oshd025_GL, snow_total_yr_oshd025_GL, \
            swe_total_Hyr_oshd025_GL = open_files('OSHD_FILES_025_new', bf, yr_in)
        # cru no lapse
        gpp_cru1km_GL, gpp0506_cru1km_GL, evap_cru1km_GL, evap_cru1km_GL_HY, melt_out_cru1km_GL, melt_doy_cru1km_GL, \
            snow_cru1km_GL, snow_yr_cru1km_GL, snow_Hyr_cru1km_GL, snow_total_cru1km_GL, snow_total_yr_cru1km_GL, \
            swe_total_Hyr_cru1km_GL = open_files('CRUJRA_FILES_noLapse_OLD', bf, yr_in)
        gpp_cru1km_HR, gpp0506_cru1km_HR, evap_cru1km_HR, evap_cru1km_HR_HY, melt_out_cru1km_HR, melt_doy_cru1km_HR, \
            snow_cru1km_HR, snow_yr_cru1km_HR, snow_Hyr_cru1km_HR, snow_total_cru1km_HR, snow_total_yr_cru1km_HR, \
            swe_total_Hyr_cru1km_HR = open_files('CRUJRA_FILES_noLapse', bf, yr_in)
        gpp_cru05_GL, gpp0506_cru05_GL, evap_cru05_GL, evap_cru05_GL_HY, melt_out_cru05_GL, melt_doy_cru05_GL, \
            snow_cru05_GL, snow_yr_cru05_GL, snow_Hyr_cru05_GL, snow_total_cru05_GL, snow_total_yr_cru05_GL, \
            swe_total_Hyr_cru05_GL = open_files('CRUJRA_FILES_05deg_cru_new', bf, yr_in)
        gpp_cru025_GL, gpp0506_cru025_GL, evap_cru025_GL, evap_cru025_GL_HY, melt_out_cru025_GL, melt_doy_cru025_GL, \
            snow_cru025_GL, snow_yr_cru025_GL, snow_Hyr_cru025_GL, snow_total_cru025_GL, snow_total_yr_cru025_GL, \
            swe_total_Hyr_cru025_GL = open_files('CRUJRA_FILES_025deg_cru_new', bf, yr_in)
        # cru plus with lapse rate
        gpp_cruL1km_GL, gpp0506_cruL1km_GL, evap_cruL1km_GL, evap_cruL1km_GL_HY, melt_out_cruL1km_GL, melt_doy_cruL1km_GL, \
            snow_cruL1km_GL, snow_yr_cruL1km_GL, snow_Hyr_cruL1km_GL, snow_total_cruL1km_GL, snow_total_yr_cruL1km_GL, \
            swe_total_Hyr_cruL1km_GL = open_files('CRUJRA_FILES_OLD', bf, yr_in)
        gpp_cruL1km_HR, gpp0506_cruL1km_HR, evap_cruL1km_HR, evap_cruL1km_HR_HY, melt_out_cruL1km_HR, melt_doy_cruL1km_HR, \
            snow_cruL1km_HR, snow_yr_cruL1km_HR, snow_Hyr_cruL1km_HR, snow_total_cruL1km_HR, snow_total_yr_cruL1km_HR, \
            swe_total_Hyr_cruL1km_HR = open_files('CRUJRA_FILES', bf, yr_in)
        gpp_cruL05_GL, gpp0506_cruL05_GL, evap_cruL05_GL, evap_cruL05_GL_HY, melt_out_cruL05_GL, melt_doy_cruL05_GL, \
            snow_cruL05_GL, snow_yr_cruL05_GL, snow_Hyr_cruL05_GL, snow_total_cruL05_GL, snow_total_yr_cruL05_GL, \
            swe_total_Hyr_cruL05_GL = open_files('CRUJRA_FILES_05deg_cru_new_lapse', bf, yr_in)
        gpp_cruL025_GL, gpp0506_cruL025_GL, evap_cruL025_GL, evap_cruL025_GL_HY, melt_out_cruL025_GL, melt_doy_cruL025_GL, \
            snow_cruL025_GL, snow_yr_cruL025_GL, snow_Hyr_cruL025_GL, snow_total_cruL025_GL, snow_total_yr_cruL025_GL, \
            swe_total_Hyr_cruL025_GL = open_files('CRUJRA_FILES_025deg_cru_new_lapse', bf, yr_in)

        proj_data = ccrs.PlateCarree()
        proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)

        fig, axs = plt.subplots(1, 2, figsize=[7, 4], subplot_kw={'projection': proj_map})
        a1 = gpp_oshd1km_HR.DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[0])
        axs[0].add_feature(cf.BORDERS, linewidth=1.2, edgecolor='darkred', alpha=1)
        axs[0].axis('off')
        axs[0].set_title('', fontweight='bold', fontsize=14)
        cbar = fig.colorbar(a1, ax=axs[0], shrink=0.6, location='bottom')
        string_units = 'GPP$_{JulAug}$ [gC m$^{-2}$ mo$^{-1}$]'
        cbar.set_label(string_units, rotation=0, labelpad=1)

        a1 = snow_oshd1km_HR.DATA.plot(transform=proj_data, cmap='viridis', add_colorbar=False, ax=axs[1])
        axs[1].add_feature(cf.BORDERS, linewidth=1.2, edgecolor='darkred', alpha=1)
        axs[1].axis('off')
        axs[1].set_title('', fontweight='bold', fontsize=14)
        cbar = fig.colorbar(a1, ax=axs[1], shrink=0.6, location='bottom')
        string_units = 'no. days w/ >2cm snow Jan-Jun'
        cbar.set_label(string_units, rotation=0, labelpad=1)
        plt.tight_layout()
        fig.savefig(bf_plots / Path('absolute_gpp_snow_ref_' + str(yr_in) + '.png'), dpi=500,
                    bbox_inches='tight', facecolor='white', transparent=False)

        # hexbin plot
        x1 = snow_oshd1km_HR.where((snow_oshd1km_HR.mask_veg == 1)).data.flatten()
        y1 = gpp_oshd1km_HR.where((gpp_oshd1km_HR.mask_veg == 1)).data.flatten()
        id1 = (~np.isnan(x1)) | (~np.isnan(y1))

        fig = plt.figure()
        ax2 = plt.subplot(111)
        hb = ax2.hexbin(x1[id1], y1[id1], gridsize=30, mincnt=10, cmap='inferno')
        cb = fig.colorbar(hb, ax=ax2, label='counts')
        ax2.set_xlabel('# Days > 2cm Snow')
        ax2.set_ylabel('GPP$_{JA}$ [gC m$^{-2}$ month$^{-1}$]')
        fig.savefig(Path(r'C:\Users\malle\Documents\paper_clm5_ch\fig_s9_hexbin.pdf'), transparent=True)

        # look at snow days b/w january and june Diffs
        string_units = 'no. days w/ >2cm snow Jan-Jun'
        string_save = 'snow_days_jan_june_'
        make_rel_plots(snow_oshd1km_HR, snow_oshd1km_GL, snow_oshd025_GL, snow_oshd05_GL, snow_cru1km_HR,
                       snow_cru1km_GL, snow_cru025_GL, snow_cru05_GL, snow_cruL1km_HR, snow_cruL1km_GL, snow_cruL025_GL,
                       snow_cruL05_GL, string_units, string_save, bf_plots)

        # look at GPP July/August Diffs
        string_units = 'GPP [gC m$^{-2}$ mo$^{-1}$]'
        string_save = 'gpp_jul_aug_'
        make_rel_plots(gpp_oshd1km_HR, gpp_oshd1km_GL, gpp_oshd025_GL, gpp_oshd05_GL, gpp_cru1km_HR, gpp_cru1km_GL,
                       gpp_cru025_GL, gpp_cru05_GL, gpp_cruL1km_HR, gpp_cruL1km_GL, gpp_cruL025_GL, gpp_cruL05_GL,
                       string_units, string_save, bf_plots)

        # look at GPP May/June Diffs
        string_save = 'gpp_may_jun_'
        make_rel_plots(gpp0506_oshd1km_HR, gpp0506_oshd1km_GL, gpp0506_oshd025_GL, gpp0506_oshd05_GL, gpp0506_cru1km_HR,
                       gpp0506_cru1km_GL, gpp0506_cru025_GL, gpp0506_cru05_GL, gpp0506_cruL1km_HR, gpp0506_cruL1km_GL,
                       gpp0506_cruL025_GL, gpp0506_cruL05_GL, string_units, string_save, bf_plots)

        # look at yearly ET Diffs
        string_units = 'ET [mm yr$^{-1}$]'
        string_save = 'ET_yearly_'
        make_rel_plots(evap_oshd1km_HR, evap_oshd1km_GL, evap_oshd025_GL, evap_oshd05_GL, evap_cru1km_HR,
                       evap_cru1km_GL, evap_cru025_GL, evap_cru05_GL, evap_cruL1km_HR, evap_cruL1km_GL, evap_cruL025_GL,
                       evap_cruL05_GL, string_units, string_save, bf_plots)

        # look at hydro-yearly ET Diffs
        string_units = 'ET [mm yr$^{-1}$]'
        string_save = 'ET_hydro_yearly_'
        make_rel_plots(evap_oshd1km_HR_HY, evap_oshd1km_GL_HY, evap_oshd025_GL_HY, evap_oshd05_GL_HY, evap_cru1km_HR_HY,
                       evap_cru1km_GL_HY, evap_cru025_GL_HY, evap_cru05_GL_HY, evap_cruL1km_HR_HY, evap_cruL1km_GL_HY,
                       evap_cruL025_GL_HY, evap_cruL05_GL_HY, string_units, string_save, bf_plots)

        # look at SWE hydro-year Diffs
        string_units = 'SWE [mm yr$^{-1}$]'
        string_save = 'SWE_hydro_yearly'
        make_rel_plots(swe_total_Hyr_oshd1km_HR, swe_total_Hyr_oshd1km_GL, swe_total_Hyr_oshd025_GL,
                       swe_total_Hyr_oshd05_GL, swe_total_Hyr_cru1km_HR, swe_total_Hyr_cru1km_GL,
                       swe_total_Hyr_cru025_GL, swe_total_Hyr_cru05_GL, swe_total_Hyr_cruL1km_HR,
                       swe_total_Hyr_cruL1km_GL, swe_total_Hyr_cruL025_GL, swe_total_Hyr_cruL05_GL, string_units,
                       string_save, bf_plots)

        # look at melt-out date
        string_units = 'Melt-out [DOY]'
        string_save = 'melt_out_doy'
        make_rel_plots(melt_doy_oshd1km_HR, melt_doy_oshd1km_GL, melt_doy_oshd025_GL, melt_doy_oshd05_GL,
                       melt_doy_cru1km_HR, melt_doy_cru1km_GL, melt_doy_cru025_GL, melt_doy_cru05_GL,
                       melt_doy_cruL1km_HR, melt_doy_cruL1km_GL, melt_doy_cruL025_GL, melt_doy_cruL05_GL, string_units,
                       string_save, bf_plots)

        # make boxplots of GPP & snow cover
        df_gpp0506 = make_df_boxplots(gpp0506_oshd1km_HR, gpp0506_oshd1km_GL, gpp0506_oshd025_GL, gpp0506_oshd05_GL,
                                      gpp0506_cru1km_HR, gpp0506_cru1km_GL, gpp0506_cru025_GL, gpp0506_cru05_GL,
                                      gpp0506_cruL1km_HR, gpp0506_cruL1km_GL, gpp0506_cruL025_GL, gpp0506_cruL05_GL)

        df_gpp0708 = make_df_boxplots(gpp_oshd1km_HR, gpp_oshd1km_GL, gpp_oshd025_GL, gpp_oshd05_GL, gpp_cru1km_HR,
                                      gpp_cru1km_GL, gpp_cru025_GL, gpp_cru05_GL, gpp_cruL1km_HR, gpp_cruL1km_GL,
                                      gpp_cruL025_GL, gpp_cruL05_GL)

        df_snow0506 = make_df_boxplots(snow_oshd1km_HR, snow_oshd1km_GL, snow_oshd025_GL, snow_oshd05_GL,
                                       snow_cru1km_HR, snow_cru1km_GL, snow_cru025_GL, snow_cru05_GL, snow_cruL1km_HR,
                                       snow_cruL1km_GL, snow_cruL025_GL, snow_cruL05_GL)

        df_melt_doy = make_df_boxplots(melt_doy_oshd1km_HR, melt_doy_oshd1km_GL, melt_doy_oshd025_GL,
                                       melt_doy_oshd05_GL, melt_doy_cru1km_HR, melt_doy_cru1km_GL, melt_doy_cru025_GL,
                                       melt_doy_cru05_GL, melt_doy_cruL1km_HR, melt_doy_cruL1km_GL, melt_doy_cruL025_GL,
                                       melt_doy_cruL05_GL)

        df_evap_hyd = make_df_boxplots(evap_oshd1km_HR_HY, evap_oshd1km_GL_HY, evap_oshd025_GL_HY, evap_oshd05_GL_HY,
                                       evap_cru1km_HR_HY, evap_cru1km_GL_HY, evap_cru025_GL_HY, evap_cru05_GL_HY,
                                       evap_cruL1km_HR_HY, evap_cruL1km_GL_HY, evap_cruL025_GL_HY, evap_cruL05_GL)

        df_swe_hyd = make_df_boxplots(swe_total_Hyr_oshd1km_HR, swe_total_Hyr_oshd1km_GL, swe_total_Hyr_oshd025_GL,
                                      swe_total_Hyr_oshd05_GL, swe_total_Hyr_cru1km_HR, swe_total_Hyr_cru1km_GL,
                                      swe_total_Hyr_cru025_GL, swe_total_Hyr_cru05_GL, swe_total_Hyr_cruL1km_HR,
                                      swe_total_Hyr_cruL1km_GL, swe_total_Hyr_cruL025_GL, swe_total_Hyr_cruL05_GL)

        sns.set_style('white', {'axes.linewidth': 0.5})
        plt.rcParams['xtick.major.size'] = 15
        plt.rcParams['xtick.major.width'] = 3
        plt.rcParams['xtick.bottom'] = True
        bf_plots_box = bf / 'snow_eco_new' / 'gpp_snow'

        fig = plt.figure()
        axes = fig.add_subplot(111)
        ax1 = sns.violinplot(data=df_gpp0708, cut=0, palette=my_pal, split=True)
        plt.setp(ax1.get_xticklabels(), rotation=20, ha='right', rotation_mode='anchor')
        axes.yaxis.grid(True)
        axes.set_ylabel('GPP$_{JulAug}$ [gC m$^{-2}$ month$^{-1}$]', fontsize=16)
        for patch in ax1.collections[2::8]:
            patch.set_alpha(0.7)
        for patch in ax1.collections[4::8]:
            patch.set_alpha(0.45)
        for patch in ax1.collections[6::8]:
            patch.set_alpha(0.2)

        for tick_label in ax1.get_xticklabels()[:4]:
            tick_label.set_color((27 / 255, 158 / 255, 119 / 255))
        for tick_label in ax1.get_xticklabels()[4:8]:
            tick_label.set_color((117 / 255, 112 / 255, 179 / 255))
        for tick_label in ax1.get_xticklabels()[8:]:
            tick_label.set_color((217 / 255, 95 / 255, 2 / 255))

        for tick_label in ax1.get_xticklabels()[1::4]:
            tick_label.set_alpha(0.7)
        for tick_label in ax1.get_xticklabels()[2::4]:
            tick_label.set_alpha(0.5)
        for tick_label in ax1.get_xticklabels()[3::4]:
            tick_label.set_alpha(0.3)
        plt.tight_layout()
        fig.savefig(bf_plots_box / Path('comp_vio_gpp0708_' + str(yr_in) + '.png'),
                    facecolor='white', transparent=False)

        fig = plt.figure(figsize=(15, 10))
        plt.suptitle(str(yr_in), fontsize=16)
        axes = fig.add_subplot(221)
        ax1 = sns.violinplot(data=df_snow0506, cut=0, palette=my_pal, split=True)
        plt.title('(a)', ha='left', x=-0.04)
        axes.tick_params(labelbottom=False)
        axes.yaxis.grid(True)
        axes.set_ylabel('no. days w/ >2cm snow Jan-Jun', fontsize=16)
        for patch in ax1.collections[2::8]:
            patch.set_alpha(0.7)
        for patch in ax1.collections[4::8]:
            patch.set_alpha(0.45)
        for patch in ax1.collections[6::8]:
            patch.set_alpha(0.2)

        axes = fig.add_subplot(222)
        ax1 = sns.violinplot(data=df_swe_hyd, cut=0, palette=my_pal, split=True)
        plt.title('(b)', ha='left', x=-0.04)
        axes.tick_params(labelbottom=False)
        axes.yaxis.grid(True)
        axes.set_ylabel('SWE$_{HydYr}$ [mm yr$^{-1}$]', fontsize=16)
        for patch in ax1.collections[2::8]:
            patch.set_alpha(0.7)
        for patch in ax1.collections[4::8]:
            patch.set_alpha(0.45)
        for patch in ax1.collections[6::8]:
            patch.set_alpha(0.2)

        axes = fig.add_subplot(223)
        ax1 = sns.violinplot(data=df_gpp0506, cut=0, palette=my_pal, split=True)
        plt.title('(c)', ha='left', x=-0.04)
        plt.setp(ax1.get_xticklabels(), rotation=20, ha='right', rotation_mode='anchor')
        axes.yaxis.grid(True)
        axes.set_ylabel('GPP$_{MayJune}$ [gC m$^{-2}$ month$^{-1}$]', fontsize=16)
        for patch in ax1.collections[2::8]:
            patch.set_alpha(0.7)
        for patch in ax1.collections[4::8]:
            patch.set_alpha(0.45)
        for patch in ax1.collections[6::8]:
            patch.set_alpha(0.2)

        for tick_label in ax1.get_xticklabels()[:4]:
            tick_label.set_color((27 / 255, 158 / 255, 119 / 255))
        for tick_label in ax1.get_xticklabels()[4:8]:
            tick_label.set_color((117 / 255, 112 / 255, 179 / 255))
        for tick_label in ax1.get_xticklabels()[8:]:
            tick_label.set_color((217 / 255, 95 / 255, 2 / 255))

        for tick_label in ax1.get_xticklabels()[1::4]:
            tick_label.set_alpha(0.7)
        for tick_label in ax1.get_xticklabels()[2::4]:
            tick_label.set_alpha(0.5)
        for tick_label in ax1.get_xticklabels()[3::4]:
            tick_label.set_alpha(0.3)

        axes = fig.add_subplot(224)
        ax1 = sns.violinplot(data=df_evap_hyd, cut=0, palette=my_pal, split=True)
        plt.title('(d)', ha='left', x=-0.04)
        plt.setp(ax1.get_xticklabels(), rotation=20, ha='right', rotation_mode='anchor')
        axes.yaxis.grid(True)
        axes.set_ylabel('ET$_{HydYr}$ [mm yr$^{-1}$]', fontsize=16)
        plt.tight_layout()
        for patch in ax1.collections[2::8]:
            patch.set_alpha(0.7)
        for patch in ax1.collections[4::8]:
            patch.set_alpha(0.45)
        for patch in ax1.collections[6::8]:
            patch.set_alpha(0.2)

        for tick_label in ax1.get_xticklabels()[:4]:
            tick_label.set_color((27 / 255, 158 / 255, 119 / 255))
        for tick_label in ax1.get_xticklabels()[4:8]:
            tick_label.set_color((117 / 255, 112 / 255, 179 / 255))
        for tick_label in ax1.get_xticklabels()[8:]:
            tick_label.set_color((217 / 255, 95 / 255, 2 / 255))

        for tick_label in ax1.get_xticklabels()[1::4]:
            tick_label.set_alpha(0.7)
        for tick_label in ax1.get_xticklabels()[2::4]:
            tick_label.set_alpha(0.5)
        for tick_label in ax1.get_xticklabels()[3::4]:
            tick_label.set_alpha(0.3)
        fig.savefig(bf_plots_box / Path('comp_vio_' + str(yr_in) + '.png'),
                    facecolor='white', transparent=False)

        # make boxplots of GPP & snow cover -> only points below 3000m -> because we set to nan based on oshd1km only need to use mask there
        df_gpp0506_3000 = make_df_boxplots(gpp0506_oshd1km_HR.where((gpp0506_oshd1km_HR.mask_3000 == 0)),
                                           gpp0506_oshd1km_GL, gpp0506_oshd025_GL, gpp0506_oshd05_GL, gpp0506_cru1km_HR,
                                           gpp0506_cru1km_GL, gpp0506_cru025_GL, gpp0506_cru05_GL, gpp0506_cruL1km_HR,
                                           gpp0506_cruL1km_GL, gpp0506_cruL025_GL, gpp0506_cruL05_GL)

        df_gpp0708_3000 = make_df_boxplots(gpp_oshd1km_HR.where((gpp0506_oshd1km_HR.mask_3000 == 0)), gpp_oshd1km_GL,
                                           gpp_oshd025_GL, gpp_oshd05_GL, gpp_cru1km_HR, gpp_cru1km_GL, gpp_cru025_GL,
                                           gpp_cru05_GL, gpp_cruL1km_HR, gpp_cruL1km_GL, gpp_cruL025_GL, gpp_cruL05_GL)

        df_snow0506_3000 = make_df_boxplots(snow_oshd1km_HR.where((gpp0506_oshd1km_HR.mask_3000 == 0)), snow_oshd1km_GL,
                                            snow_oshd025_GL, snow_oshd05_GL, snow_cru1km_HR, snow_cru1km_GL,
                                            snow_cru025_GL, snow_cru05_GL, snow_cruL1km_HR, snow_cruL1km_GL,
                                            snow_cruL025_GL, snow_cruL05_GL)

        df_melt_doy_3000 = make_df_boxplots(melt_doy_oshd1km_HR.where((gpp0506_oshd1km_HR.mask_3000 == 0)),
                                            melt_doy_oshd1km_GL, melt_doy_oshd025_GL, melt_doy_oshd05_GL,
                                            melt_doy_cru1km_HR, melt_doy_cru1km_GL, melt_doy_cru025_GL,
                                            melt_doy_cru05_GL, melt_doy_cruL1km_HR, melt_doy_cruL1km_GL,
                                            melt_doy_cruL025_GL, melt_doy_cruL05_GL)

        df_evap_hyd_3000 = make_df_boxplots(evap_oshd1km_HR_HY.where((gpp0506_oshd1km_HR.mask_3000 == 0)),
                                            evap_oshd1km_GL_HY, evap_oshd025_GL_HY, evap_oshd05_GL_HY,
                                            evap_cru1km_HR_HY, evap_cru1km_GL_HY, evap_cru025_GL_HY, evap_cru05_GL_HY,
                                            evap_cruL1km_HR_HY, evap_cruL1km_GL_HY, evap_cruL025_GL_HY, evap_cruL05_GL)

        df_swe_hyd_3000 = make_df_boxplots(swe_total_Hyr_oshd1km_HR.where((gpp0506_oshd1km_HR.mask_3000 == 0)),
                                           swe_total_Hyr_oshd1km_GL, swe_total_Hyr_oshd025_GL, swe_total_Hyr_oshd05_GL,
                                           swe_total_Hyr_cru1km_HR, swe_total_Hyr_cru1km_GL, swe_total_Hyr_cru025_GL,
                                           swe_total_Hyr_cru05_GL, swe_total_Hyr_cruL1km_HR, swe_total_Hyr_cruL1km_GL,
                                           swe_total_Hyr_cruL025_GL, swe_total_Hyr_cruL05_GL)

        fig = plt.figure(figsize=(22, 11))
        plt.suptitle(str(yr_in) + ', <3000m', fontsize=16)
        axes = fig.add_subplot(231)
        ax1 = sns.violinplot(data=df_melt_doy_3000, cut=0, palette=my_pal, split=True)
        axes.tick_params(labelbottom=False)
        axes.yaxis.grid(True)
        axes.set_ylabel('Melt-out [DOY]', fontsize=16)
        for patch in ax1.collections[2::8]:
            patch.set_alpha(0.7)
        for patch in ax1.collections[4::8]:
            patch.set_alpha(0.45)
        for patch in ax1.collections[6::8]:
            patch.set_alpha(0.2)

        axes = fig.add_subplot(232)
        ax1 = sns.violinplot(data=df_snow0506_3000, cut=0, palette=my_pal, split=True)
        axes.tick_params(labelbottom=False)
        axes.yaxis.grid(True)
        axes.set_ylabel('no. days w/ >2cm snow Jan-Jun', fontsize=16)
        for patch in ax1.collections[2::8]:
            patch.set_alpha(0.7)
        for patch in ax1.collections[4::8]:
            patch.set_alpha(0.45)
        for patch in ax1.collections[6::8]:
            patch.set_alpha(0.2)

        axes = fig.add_subplot(233)
        ax1 = sns.violinplot(data=df_swe_hyd_3000, cut=0, palette=my_pal, split=True)
        axes.tick_params(labelbottom=False)
        axes.yaxis.grid(True)
        axes.set_ylabel('SWE$_{HydYr}$ [mm yr$^{-1}$]', fontsize=16)
        for patch in ax1.collections[2::8]:
            patch.set_alpha(0.7)
        for patch in ax1.collections[4::8]:
            patch.set_alpha(0.45)
        for patch in ax1.collections[6::8]:
            patch.set_alpha(0.2)

        axes = fig.add_subplot(234)
        ax1 = sns.violinplot(data=df_gpp0506_3000, cut=0, palette=my_pal, split=True)
        plt.setp(ax1.get_xticklabels(), rotation=20, ha='right', rotation_mode='anchor')
        axes.yaxis.grid(True)
        axes.set_ylabel('GPP$_{MayJune}$ [gC m$^{-2}$ month$^{-1}$]', fontsize=16)
        for patch in ax1.collections[2::8]:
            patch.set_alpha(0.7)
        for patch in ax1.collections[4::8]:
            patch.set_alpha(0.45)
        for patch in ax1.collections[6::8]:
            patch.set_alpha(0.2)

        for tick_label in ax1.get_xticklabels()[:4]:
            tick_label.set_color((27 / 255, 158 / 255, 119 / 255))
        for tick_label in ax1.get_xticklabels()[4:8]:
            tick_label.set_color((117 / 255, 112 / 255, 179 / 255))
        for tick_label in ax1.get_xticklabels()[8:]:
            tick_label.set_color((217 / 255, 95 / 255, 2 / 255))

        for tick_label in ax1.get_xticklabels()[1::4]:
            tick_label.set_alpha(0.7)
        for tick_label in ax1.get_xticklabels()[2::4]:
            tick_label.set_alpha(0.5)
        for tick_label in ax1.get_xticklabels()[3::4]:
            tick_label.set_alpha(0.3)

        axes = fig.add_subplot(235)
        ax1 = sns.violinplot(data=df_gpp0708_3000, cut=0, palette=my_pal, split=True)
        plt.setp(ax1.get_xticklabels(), rotation=20, ha='right', rotation_mode='anchor')
        axes.yaxis.grid(True)
        axes.set_ylabel('GPP$_{JulAug}$ [gC m$^{-2}$ month$^{-1}$]', fontsize=16)
        for patch in ax1.collections[2::8]:
            patch.set_alpha(0.7)
        for patch in ax1.collections[4::8]:
            patch.set_alpha(0.45)
        for patch in ax1.collections[6::8]:
            patch.set_alpha(0.2)

        for tick_label in ax1.get_xticklabels()[:4]:
            tick_label.set_color((27 / 255, 158 / 255, 119 / 255))
        for tick_label in ax1.get_xticklabels()[4:8]:
            tick_label.set_color((117 / 255, 112 / 255, 179 / 255))
        for tick_label in ax1.get_xticklabels()[8:]:
            tick_label.set_color((217 / 255, 95 / 255, 2 / 255))
        for tick_label in ax1.get_xticklabels()[1::4]:
            tick_label.set_alpha(0.7)
        for tick_label in ax1.get_xticklabels()[2::4]:
            tick_label.set_alpha(0.5)
        for tick_label in ax1.get_xticklabels()[3::4]:
            tick_label.set_alpha(0.3)

        axes = fig.add_subplot(236)
        ax1 = sns.violinplot(data=df_evap_hyd_3000, cut=0, palette=my_pal, split=True)
        plt.setp(ax1.get_xticklabels(), rotation=20, ha='right', rotation_mode='anchor')
        axes.yaxis.grid(True)
        axes.set_ylabel('ET$_{HydYr}$ [mm yr$^{-1}$]', fontsize=16)
        plt.tight_layout()
        for patch in ax1.collections[2::8]:
            patch.set_alpha(0.7)
        for patch in ax1.collections[4::8]:
            patch.set_alpha(0.45)
        for patch in ax1.collections[6::8]:
            patch.set_alpha(0.2)

        for tick_label in ax1.get_xticklabels()[:4]:
            tick_label.set_color((27 / 255, 158 / 255, 119 / 255))
        for tick_label in ax1.get_xticklabels()[4:8]:
            tick_label.set_color((117 / 255, 112 / 255, 179 / 255))
        for tick_label in ax1.get_xticklabels()[8:]:
            tick_label.set_color((217 / 255, 95 / 255, 2 / 255))
        for tick_label in ax1.get_xticklabels()[1::4]:
            tick_label.set_alpha(0.7)
        for tick_label in ax1.get_xticklabels()[2::4]:
            tick_label.set_alpha(0.5)
        for tick_label in ax1.get_xticklabels()[3::4]:
            tick_label.set_alpha(0.3)

        fig.savefig(bf_plots_box / Path('comp_vio_below3000_' + str(yr_in) + '.png'),
                    facecolor='white', transparent=False)

