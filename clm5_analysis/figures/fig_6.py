# -*- coding: utf-8 -*-
"""
Desc: Script to create difference plots as well as correlation shown in Figure 6.
Created on 13.07.23 15:34
@author: malle
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib
import pandas as pd
import seaborn as sns
import scipy
from sklearn.metrics import mean_squared_error
import platform


def calc_stats(data, x, y):
    _, _, rvalue, pvalue, _ = scipy.stats.linregress(x=data[x], y=data[y])
    rmse = mean_squared_error(data[x], data[y], squared=False)
    return rmse, pvalue, rvalue ** 2


def annotate(ax, data, x, y, col_in, y_loc_in):
    _, _, rvalue, _, _ = scipy.stats.linregress(x=data[x], y=data[y])
    rmse = mean_squared_error(data[x], data[y], squared=False)
    ax.text(.02, y_loc_in, f'R$^{2}$={rvalue ** 2:.2f}, RMSE={rmse:.2f}', transform=ax.transAxes, color=col_in)


def open_files(run_in, bf_in, year_in):
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

    gpp = xr.open_dataset(bf_in / run_in / file_in_gpp)
    gpp0506 = xr.open_dataset(bf_in / run_in / file_in_gpp_0506)
    snow = xr.open_dataset(bf_in / run_in / file_in_snow)
    snow_yr = xr.open_dataset(bf_in / run_in / file_in_snow_yr)
    snow_hyr = xr.open_dataset(bf_in / run_in / file_in_snow_hydroyr)
    evap = xr.open_dataset(bf_in / run_in / file_in_evap_yr)
    evap_hyr = xr.open_dataset(bf_in / run_in / file_in_evap_hydroyr)
    snow_total = xr.open_dataset(bf_in / run_in / file_in_snow_total)
    snow_total_yr = xr.open_dataset(bf_in / run_in / file_in_snow_total_yr)
    swe_total_hyr = xr.open_dataset(bf_in / run_in / file_in_swe_total)
    swe_total_hyr = swe_total_hyr.where(snow_hyr.DATA > 50, np.nan)

    melt_out = xr.open_dataset(bf_in / run_in / file_in_melt_date)
    melt_doy = xr.open_dataset(bf_in / run_in / file_in_melt_doy)
    # set melt_out dates to nan when < 50 snow days per hydrological year or when later than July 1st
    melt_out = melt_out.where(snow_hyr.DATA > 50, np.datetime64('NaT'))
    melt_doy = melt_doy.where(snow_hyr.DATA > 50, np.nan)
    melt_doy = melt_doy.where(melt_doy.DATA < 183, np.nan)
    melt_doy = melt_doy.where(melt_doy.DATA > 31, np.nan)

    return (gpp, gpp0506, evap, evap_hyr, melt_out, melt_doy,
            snow, snow_yr, snow_hyr, snow_total, snow_total_yr, swe_total_hyr)


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


if __name__ == '__main__':

    if platform.system() == 'Windows':
        bf = Path('L:\malle\CLM5_CH')
    else:
        bf = Path('/home/lud11/malle/CLM5_CH')

    bf_plots = bf / 'snow_eco_new'

    proj_data = ccrs.PlateCarree()
    proj_map = ccrs.UTM(zone=32, southern_hemisphere=False)

    comp_all_text = ['OSHD1kmHR-OSHD1kmGL', 'OSHD1kmHR-CRU1kmHR', 'OSHD1kmHR-CRU1kmGL']
    my_pal = {'< 1000m': 'yellow', '1000 - 2000m': 'blue', '> 2000m': 'grey'}
    values = ['< 1000m', '1000 - 2000m', '> 2000m']

    yr_all = [2016, 2017, 2018, 2019]
    for yr_in in yr_all:
        print(yr_in)
        (gpp_oshd1km_HR, gpp0506_oshd1km_HR, evap_oshd1km_HR, evap_oshd1km_HR_HY, melt_out_oshd1km_HR,
         melt_doy_oshd1km_HR, snow_oshd1km_HR, snow_yr_oshd1km_HR, snow_Hyr_oshd1km_HR, snow_total_oshd1km_HR,
         snow_total_yr_oshd1km_HR, swe_total_Hyr_oshd1km_HR) = open_files('OSHD_FILES', bf, yr_in)

        (gpp_oshd1km_GL, gpp0506_oshd1km_GL, evap_oshd1km_GL, evap_oshd1km_GL_HY, melt_out_oshd1km_GL,
         melt_doy_oshd1km_GL, snow_oshd1km_GL, snow_yr_oshd1km_GL, snow_Hyr_oshd1km_GL, snow_total_oshd1km_GL,
         snow_total_yr_oshd1km_GL, swe_total_Hyr_oshd1km_GL) = open_files('OSHD_FILES_OLD', bf, yr_in)

        gpp_cru1km_GL, gpp0506_cru1km_GL, evap_cru1km_GL, evap_cru1km_GL_HY, melt_out_cru1km_GL, melt_doy_cru1km_GL, \
            snow_cru1km_GL, snow_yr_cru1km_GL, snow_Hyr_cru1km_GL, snow_total_cru1km_GL, snow_total_yr_cru1km_GL, \
            swe_total_Hyr_cru1km_GL = open_files('CRUJRA_FILES_noLapse_OLD', bf, yr_in)
        gpp_cru1km_HR, gpp0506_cru1km_HR, evap_cru1km_HR, evap_cru1km_HR_HY, melt_out_cru1km_HR, melt_doy_cru1km_HR, \
            snow_cru1km_HR, snow_yr_cru1km_HR, snow_Hyr_cru1km_HR, snow_total_cru1km_HR, snow_total_yr_cru1km_HR, \
            swe_total_Hyr_cru1km_HR = open_files('CRUJRA_FILES_noLapse', bf, yr_in)

        for comp_text in comp_all_text:
            # first load everything
            if comp_text == 'OSHD1kmHR-OSHD1kmGL':
                bf_plots_in = bf_plots / 'oshd_effect_surfdata'
                melt_1 = melt_out_oshd1km_HR.where((melt_out_oshd1km_HR.mask_veg == 1)).DATA
                melt_2 = melt_out_oshd1km_GL.where((melt_out_oshd1km_GL.mask_veg == 1)).DATA
                gpp_1 = gpp_oshd1km_HR.where((gpp_oshd1km_HR.mask_veg == 1)).DATA
                gpp_2 = gpp_oshd1km_GL.where((gpp_oshd1km_GL.mask_veg == 1)).DATA
                gpp_1_0506 = gpp0506_oshd1km_HR.where((gpp0506_oshd1km_HR.mask_veg == 1)).DATA
                gpp_2_0506 = gpp0506_oshd1km_GL.where((gpp0506_oshd1km_GL.mask_veg == 1)).DATA
                evapo_1 = evap_oshd1km_HR_HY.where((evap_oshd1km_HR_HY.mask_veg == 1)).DATA
                evapo_2 = evap_oshd1km_GL_HY.where((evap_oshd1km_GL_HY.mask_veg == 1)).DATA
                swe_total_1 = swe_total_Hyr_oshd1km_HR.where((swe_total_Hyr_oshd1km_HR.mask_veg == 1)).DATA
                swe_total_2 = swe_total_Hyr_oshd1km_GL.where((swe_total_Hyr_oshd1km_GL.mask_veg == 1)).DATA
                snow_days_1 = snow_oshd1km_HR.where((snow_oshd1km_HR.mask_veg == 1)).DATA
                snow_days_2 = snow_oshd1km_GL.where((snow_oshd1km_GL.mask_veg == 1)).DATA
                melt_doy_1 = melt_doy_oshd1km_HR.where((melt_doy_oshd1km_HR.mask_veg == 1)).DATA
                melt_doy_2 = melt_doy_oshd1km_GL.where((melt_doy_oshd1km_GL.mask_veg == 1)).DATA
            elif comp_text == 'OSHD1kmHR-CRU1kmHR':
                bf_plots_in = bf_plots / 'oshd_cru_effect_meteo'
                melt_1 = melt_out_oshd1km_HR.where((melt_out_oshd1km_HR.mask_veg == 1)).DATA
                melt_2 = melt_out_cru1km_HR.where((melt_out_oshd1km_GL.mask_veg == 1)).DATA
                gpp_1 = gpp_oshd1km_HR.where((gpp_oshd1km_HR.mask_veg == 1)).DATA
                gpp_2 = gpp_cru1km_HR.where((gpp_oshd1km_GL.mask_veg == 1)).DATA
                gpp_1_0506 = gpp0506_oshd1km_HR.where((gpp0506_oshd1km_HR.mask_veg == 1)).DATA
                gpp_2_0506 = gpp0506_cru1km_HR.where((gpp0506_oshd1km_GL.mask_veg == 1)).DATA
                evapo_1 = evap_oshd1km_HR_HY.where((evap_oshd1km_HR_HY.mask_veg == 1)).DATA
                evapo_2 = evap_cru1km_HR.where((evap_oshd1km_GL_HY.mask_veg == 1)).DATA
                swe_total_1 = swe_total_Hyr_oshd1km_HR.where((swe_total_Hyr_oshd1km_HR.mask_veg == 1)).DATA
                swe_total_2 = swe_total_Hyr_cru1km_HR.where((swe_total_Hyr_oshd1km_GL.mask_veg == 1)).DATA
                snow_days_1 = snow_oshd1km_HR.where((snow_oshd1km_HR.mask_veg == 1)).DATA
                snow_days_2 = snow_cru1km_HR.where((snow_oshd1km_GL.mask_veg == 1)).DATA
                melt_doy_1 = melt_doy_oshd1km_HR.where((melt_doy_oshd1km_HR.mask_veg == 1)).DATA
                melt_doy_2 = melt_doy_cru1km_HR.where((melt_doy_cru1km_HR.mask_veg == 1)).DATA
            elif comp_text == 'OSHD1kmHR-CRU1kmGL':
                bf_plots_in = bf_plots / 'oshd_cru_effect_meteo_surfdata'
                melt_1 = melt_out_oshd1km_HR.where((melt_out_oshd1km_HR.mask_veg == 1)).DATA
                melt_2 = melt_out_cru1km_GL.where((melt_out_oshd1km_GL.mask_veg == 1)).DATA
                gpp_1 = gpp_oshd1km_HR.where((gpp_oshd1km_HR.mask_veg == 1)).DATA
                gpp_2 = gpp_cru1km_GL.where((gpp_oshd1km_GL.mask_veg == 1)).DATA
                gpp_1_0506 = gpp0506_oshd1km_HR.where((gpp0506_oshd1km_HR.mask_veg == 1)).DATA
                gpp_2_0506 = gpp0506_cru1km_GL.where((gpp0506_oshd1km_GL.mask_veg == 1)).DATA
                evapo_1 = evap_oshd1km_HR_HY.where((evap_oshd1km_HR_HY.mask_veg == 1)).DATA
                evapo_2 = evap_cru1km_GL.where((evap_oshd1km_GL_HY.mask_veg == 1)).DATA
                swe_total_1 = swe_total_Hyr_oshd1km_HR.where((swe_total_Hyr_oshd1km_HR.mask_veg == 1)).DATA
                swe_total_2 = swe_total_Hyr_cru1km_GL.where((swe_total_Hyr_oshd1km_GL.mask_veg == 1)).DATA
                snow_days_1 = snow_oshd1km_HR.where((snow_oshd1km_HR.mask_veg == 1)).DATA
                snow_days_2 = snow_cru1km_GL.where((snow_oshd1km_GL.mask_veg == 1)).DATA
                melt_doy_1 = melt_doy_oshd1km_HR.where((melt_doy_oshd1km_HR.mask_veg == 1)).DATA
                melt_doy_2 = melt_doy_cru1km_GL.where((melt_doy_cru1km_GL.mask_veg == 1)).DATA
            else:
                sys.exit('comp_text must be OSHD1kmHR-OSHD1kmGL or OSHD1kmHR-CRU1kmHR')

            # difference in melt-out
            diff_melt = (melt_1.data.flatten() -
                         melt_2.data.flatten()).astype('timedelta64[D]') / np.timedelta64(1, 'D')
            diff_melt_spatial = (melt_1 - melt_2)
            a1 = diff_melt_spatial.data.astype('timedelta64[D]') / np.timedelta64(1, 'D')
            a1[a1 < -365] = np.nan
            a1[a1 > 365] = np.nan
            diff_melt[diff_melt < -365] = np.nan
            diff_melt[diff_melt > 365] = np.nan
            diff_melt_spatial.data = a1

        # difference in GPP (july/august)
            diff_gpp = gpp_2.data.flatten() - gpp_1.data.flatten()
            diff_gpp_spatial = gpp_2 - gpp_1
            diff_gpp_1000 = (gpp_2.where((gpp_2.mask_1000 == 1)).data.flatten() -
                             gpp_1.where((gpp_1.mask_1000 == 1)).data.flatten())
            diff_gpp_2000 = (gpp_2.where((gpp_2.mask_2000 == 1)).data.flatten() -
                             gpp_1.where((gpp_1.mask_2000 == 1)).data.flatten())
            diff_gpp_3000 = (gpp_2.where((gpp_2.mask_3000 == 1)).data.flatten() -
                             gpp_1.where((gpp_1.mask_3000 == 1)).data.flatten())

            # difference in GPP (may/june)
            diff_gpp_0506 = gpp_2_0506.data.flatten() - gpp_1_0506.data.flatten()
            diff_gpp_spatial_0506 = gpp_2_0506 - gpp_1_0506
            diff_gpp_0506_1000 = (gpp_2_0506.where((gpp_2.mask_1000 == 1)).data.flatten() -
                                  gpp_1_0506.where((gpp_1_0506.mask_1000 == 1)).data.flatten())
            diff_gpp_0506_2000 = (gpp_2_0506.where((gpp_2.mask_2000 == 1)).data.flatten() -
                                  gpp_1_0506.where((gpp_1_0506.mask_2000 == 1)).data.flatten())
            diff_gpp_0506_3000 = (gpp_2_0506.where((gpp_2.mask_3000 == 1)).data.flatten() -
                                  gpp_1_0506.where((gpp_1_0506.mask_3000 == 1)).data.flatten())

            # difference in evapo (hydro)
            diff_evapo = evapo_2.data.flatten() - evapo_1.data.flatten()
            diff_evapo_spatial = evapo_2 - evapo_1
            diff_evapo_1000 = (evapo_2.where((evapo_2.mask_1000 == 1)).data.flatten() -
                               evapo_1.where((evapo_1.mask_1000 == 1)).data.flatten())
            diff_evapo_2000 = (evapo_2.where((evapo_2.mask_2000 == 1)).data.flatten() -
                               evapo_1.where((evapo_1.mask_2000 == 1)).data.flatten())
            diff_evapo_3000 = (evapo_2.where((evapo_2.mask_3000 == 1)).data.flatten() -
                               evapo_1.where((evapo_1.mask_3000 == 1)).data.flatten())

            # difference in snow total (hydro)
            diff_swe_total = swe_total_2.data.flatten() - swe_total_1.data.flatten()
            diff_swe_total_spatial = swe_total_2 - swe_total_1
            diff_swe_total_1000 = (swe_total_2.where((swe_total_2.mask_1000 == 1)).data.flatten() -
                                   swe_total_1.where((swe_total_1.mask_1000 == 1)).data.flatten())
            diff_swe_total_2000 = (swe_total_2.where((swe_total_2.mask_2000 == 1)).data.flatten() -
                                   swe_total_1.where((swe_total_1.mask_2000 == 1)).data.flatten())
            diff_swe_total_3000 = (swe_total_2.where((swe_total_2.mask_3000 == 1)).data.flatten() -
                                   swe_total_1.where((swe_total_1.mask_3000 == 1)).data.flatten())

            # difference in snow days
            diff_snow_days = snow_days_2.data.flatten() - snow_days_1.data.flatten()
            diff_snow_days_spatial = snow_days_2 - snow_days_1
            diff_snow_days_1000 = (snow_days_2.where((snow_days_2.mask_1000 == 1)).data.flatten() -
                                   snow_days_1.where((snow_days_1.mask_1000 == 1)).data.flatten())
            diff_snow_days_2000 = (snow_days_2.where((snow_days_2.mask_2000 == 1)).data.flatten() -
                                   snow_days_1.where((snow_days_1.mask_2000 == 1)).data.flatten())
            diff_snow_days_3000 = (snow_days_2.where((snow_days_2.mask_3000 == 1)).data.flatten() -
                                   snow_days_1.where((snow_days_1.mask_3000 == 1)).data.flatten())

            # difference in melt-out days
            diff_melt_doy = melt_doy_2.data.flatten() - melt_doy_1.data.flatten()
            diff_melt_doy_spatial = melt_doy_2 - melt_doy_1
            diff_melt_doy_1000 = (melt_doy_2.where((melt_doy_2.mask_1000 == 1)).data.flatten() -
                                  melt_doy_1.where((melt_doy_1.mask_1000 == 1)).data.flatten())
            diff_melt_doy_2000 = (melt_doy_2.where((melt_doy_2.mask_2000 == 1)).data.flatten() -
                                  melt_doy_1.where((melt_doy_1.mask_2000 == 1)).data.flatten())
            diff_melt_doy_3000 = (melt_doy_2.where((melt_doy_2.mask_3000 == 1)).data.flatten() -
                                  melt_doy_1.where((melt_doy_1.mask_3000 == 1)).data.flatten())

            df_diff_melt_doy = pd.Series(diff_melt_doy, name='diff_melt_doy')
            df_diff_snow_days = pd.Series(diff_snow_days, name='diff_snow_days')
            df_diff_swe_total = pd.Series(diff_swe_total, name='diff_swe_total')
            df_diff_evapo = pd.Series(diff_evapo, name='diff_evapo')
            df_1000 = pd.Series(melt_doy_1.mask_1000.data.flatten(), name='less1000m')
            df_2000 = pd.Series(melt_doy_1.mask_2000.data.flatten(), name='less2000m')
            df_3000 = pd.Series(melt_doy_1.mask_3000.data.flatten(), name='more2000m')
            df_diff_gpp = pd.Series(diff_gpp_0506, name=r'diff_gpp0506')
            df_diff_gpp0708 = pd.Series(diff_gpp, name=r'diff_gpp0708')

            df_all = pd.concat([df_diff_melt_doy, df_1000, df_2000, df_3000, df_diff_gpp, df_diff_gpp0708,
                                df_diff_evapo, df_diff_swe_total, df_diff_snow_days], axis=1)

            conditions = [(df_all.less1000m == 1), (df_all.less2000m == 1), (df_all.more2000m == 1)]
            df_all['Elevation'] = np.select(conditions, values)
            df_all_evap = df_all.copy()  # I don't want nans here
            df_all_evap.loc[(np.isnan(df_all['diff_evapo']))] = np.nan
            df_all.loc[(np.isnan(df_all['diff_melt_doy']))] = np.nan

            # melt-out days vs gpp 0506
            joint = sns.jointplot(data=df_all, x="diff_melt_doy", y="diff_gpp0506", hue="Elevation", palette=my_pal,
                                  edgecolor='dimgrey', joint_kws=dict(alpha=0.6))
            ax_1000 = sns.regplot(data=df_all[df_all.Elevation == '< 1000m'], x="diff_melt_doy", y="diff_gpp0506",
                                  scatter=False, color='#CDCD00', ax=joint.ax_joint)
            ax_2000 = sns.regplot(data=df_all[df_all.Elevation == '1000 - 2000m'], x="diff_melt_doy", y="diff_gpp0506",
                                  scatter=False, color='darkblue', ax=joint.ax_joint)
            ax_3000 = sns.regplot(data=df_all[df_all.Elevation == '> 2000m'], x="diff_melt_doy", y="diff_gpp0506",
                                  scatter=False, color='dimgray', ax=joint.ax_joint)
            ax_all = sns.regplot(data=df_all, x="diff_melt_doy", y="diff_gpp0506", scatter=False,
                                 color='red', line_kws={'linewidth': 2.5, 'linestyle': '--'}, ax=joint.ax_joint)
            joint.set_axis_labels('$\Delta$ Melt-out [Days]', '$\Delta$ GPP$_{MayJune}$ [gC m-2]')
            annotate(ax_1000, data=df_all[df_all.Elevation == '< 1000m'], x='diff_melt_doy', y='diff_gpp0506',
                     col_in='#CDCD00', y_loc_in=1)
            annotate(ax_2000, data=df_all[df_all.Elevation == '1000 - 2000m'], x='diff_melt_doy', y='diff_gpp0506',
                     col_in='darkblue', y_loc_in=0.95)
            annotate(ax_3000, data=df_all[df_all.Elevation == '> 2000m'], x='diff_melt_doy', y='diff_gpp0506',
                     col_in='dimgrey', y_loc_in=0.9)
            annotate(ax_all, data=df_all.dropna(), x='diff_melt_doy', y='diff_gpp0506',
                     col_in='red', y_loc_in=0.85)
            joint.fig.tight_layout()
            plt.savefig(bf_plots_in / Path("scatter_hist_melt_gpp0506_"+str(yr_in)+'.png'), dpi=500)
            plt.close()

            # snow days vs gpp 0506
            joint = sns.jointplot(data=df_all, x="diff_snow_days", y="diff_gpp0506", hue="Elevation", palette=my_pal,
                                  edgecolor='dimgrey', joint_kws=dict(alpha=0.6), height=4)
            ax_1000 = sns.regplot(data=df_all[df_all.Elevation == '< 1000m'], x="diff_snow_days", y="diff_gpp0506",
                                  scatter=False, color='#CDCD00', ax=joint.ax_joint)
            ax_2000 = sns.regplot(data=df_all[df_all.Elevation == '1000 - 2000m'], x="diff_snow_days", y="diff_gpp0506",
                                  scatter=False, color='darkblue', ax=joint.ax_joint)
            ax_3000 = sns.regplot(data=df_all[df_all.Elevation == '> 2000m'], x="diff_snow_days", y="diff_gpp0506",
                                  scatter=False, color='dimgray', ax=joint.ax_joint)
            ax_all = sns.regplot(data=df_all, x="diff_snow_days", y="diff_gpp0506", scatter=False, color='red',
                                 line_kws={'linewidth': 2.5, 'linestyle': '--'}, ax=joint.ax_joint)
            joint.set_axis_labels('$\Delta$ #snow-days > 2cm', '$\Delta$ GPP$_{MayJune}$ [gC m-2]')
            annotate(ax_1000, data=df_all[df_all.Elevation == '< 1000m'], x='diff_snow_days', y='diff_gpp0506',
                     col_in='#CDCD00', y_loc_in=0.16)
            annotate(ax_2000, data=df_all[df_all.Elevation == '1000 - 2000m'], x='diff_snow_days', y='diff_gpp0506',
                     col_in='darkblue', y_loc_in=0.11)
            annotate(ax_3000, data=df_all[df_all.Elevation == '> 2000m'], x='diff_snow_days', y='diff_gpp0506',
                     col_in='dimgrey', y_loc_in=0.06)
            annotate(ax_all, data=df_all.dropna(), x='diff_snow_days', y='diff_gpp0506',
                     col_in='red', y_loc_in=0.01)
            sns.move_legend(joint.ax_joint, title=None, ncol=1, markerscale=0.5,  fontsize=7, loc="upper right")
            joint.fig.tight_layout()
            plt.savefig(bf_plots_in / Path("scatter_hist_snow_gpp0506_"+str(yr_in)+'_v1.png'), dpi=500)
            plt.close()

            # snow days vs gpp 0708
            joint = sns.jointplot(data=df_all, x="diff_snow_days", y="diff_gpp0708", hue="Elevation", palette=my_pal,
                                  edgecolor='dimgrey', joint_kws=dict(alpha=0.6))
            ax_1000 = sns.regplot(data=df_all[df_all.Elevation == '< 1000m'], x="diff_snow_days", y="diff_gpp0708",
                                  scatter=False, color='#CDCD00', ax=joint.ax_joint)
            ax_2000 = sns.regplot(data=df_all[df_all.Elevation == '1000 - 2000m'], x="diff_snow_days", y="diff_gpp0708",
                                  scatter=False, color='darkblue', ax=joint.ax_joint)
            ax_3000 = sns.regplot(data=df_all[df_all.Elevation == '> 2000m'], x="diff_snow_days", y="diff_gpp0708",
                                  scatter=False, color='dimgray', ax=joint.ax_joint)
            ax_all = sns.regplot(data=df_all, x="diff_snow_days", y="diff_gpp0708", scatter=False, color='red',
                                 line_kws={'linewidth': 2.5, 'linestyle': '--'}, ax=joint.ax_joint)
            joint.set_axis_labels('$\Delta$ #snow-days > 2cm', '$\Delta$ GPP$_{JulyAugust}$ [gC m-2]')
            annotate(ax_1000, data=df_all[df_all.Elevation == '< 1000m'], x='diff_snow_days', y='diff_gpp0708',
                     col_in='#CDCD00', y_loc_in=1)
            annotate(ax_2000, data=df_all[df_all.Elevation == '1000 - 2000m'], x='diff_snow_days', y='diff_gpp0708',
                     col_in='darkblue', y_loc_in=0.95)
            annotate(ax_3000, data=df_all[df_all.Elevation == '> 2000m'], x='diff_snow_days', y='diff_gpp0708',
                     col_in='dimgrey', y_loc_in=0.9)
            annotate(ax_all, data=df_all.dropna(), x='diff_snow_days', y='diff_gpp0708',
                     col_in='red', y_loc_in=0.85)
            joint.fig.tight_layout()
            plt.savefig(bf_plots_in / Path("scatter_hist_snow_gpp0708_"+str(yr_in)+'.png'), dpi=500)
            plt.close()

            # evap vs snow total
            joint = sns.jointplot(data=df_all_evap, x="diff_swe_total", y="diff_evapo", hue="Elevation", palette=my_pal,
                                  edgecolor='dimgrey', joint_kws=dict(alpha=0.6), height=4)
            ax_1000 = sns.regplot(data=df_all_evap[df_all_evap.Elevation == '< 1000m'], x="diff_swe_total",
                                  y="diff_evapo", scatter=False, color='#CDCD00', ax=joint.ax_joint)
            ax_2000 = sns.regplot(data=df_all_evap[df_all_evap.Elevation == '1000 - 2000m'], x="diff_swe_total",
                                  y="diff_evapo", scatter=False, color='darkblue', ax=joint.ax_joint)
            ax_3000 = sns.regplot(data=df_all_evap[df_all_evap.Elevation == '> 2000m'], x="diff_swe_total",
                                  y="diff_evapo", scatter=False, color='dimgray', ax=joint.ax_joint)
            ax_all = sns.regplot(data=df_all_evap, x="diff_swe_total", y="diff_evapo", scatter=False, color='red',
                                 line_kws={'linewidth': 2.5, 'linestyle': '--'}, ax=joint.ax_joint)
            joint.set_axis_labels('$\Delta$ SWE [mm yr-1]', '$\Delta$ ET [mm yr-1]')
            annotate(ax_1000, data=df_all[df_all.Elevation == '< 1000m'].dropna(), x='diff_swe_total', y='diff_evapo',
                     col_in='#CDCD00', y_loc_in=0.16)
            annotate(ax_2000, data=df_all[df_all.Elevation == '1000 - 2000m'].dropna(), x='diff_swe_total',
                     y='diff_evapo', col_in='darkblue', y_loc_in=0.11)
            annotate(ax_3000, data=df_all[df_all.Elevation == '> 2000m'].dropna(), x='diff_swe_total', y='diff_evapo',
                     col_in='dimgrey', y_loc_in=0.06)
            annotate(ax_all, data=df_all.dropna(), x='diff_swe_total', y='diff_evapo',
                     col_in='red', y_loc_in=0.01)
            sns.move_legend(joint.ax_joint, title=None, ncol=1, markerscale=0.5,  fontsize=7, loc="upper right")
            joint.fig.tight_layout()
            plt.savefig(bf_plots_in / Path("scatter_hist_swe_et_"+str(yr_in)+'_v3.png'), dpi=300)
            plt.close()

            # spatial comparison
            fig1 = plt.figure(figsize=[9, 6])  #
            ax0 = plt.subplot(221, projection=proj_map)
            diff_gpp_spatial_0506.plot(transform=proj_data, cmap='RdBu_r', add_colorbar=True,
                                       subplot_kws={'projection': proj_map},
                                       cbar_kwargs={"label": "$\Delta$ GPP$_{MayJune}$ [gC m-2]", "shrink": 0.6})
            ax0.add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
            plt.title('')
            plt.axis('off')
            ax1 = plt.subplot(222, projection=proj_map)
            diff_snow_days_spatial.plot(transform=proj_data, cmap='RdBu_r', add_colorbar=True,
                                        subplot_kws={'projection': proj_map},
                                        cbar_kwargs={"label": "$\Delta$ #snow-days > 2cm", "shrink": 0.6})
            ax1.add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
            plt.title('')
            plt.axis('off')

            ax0 = plt.subplot(223, projection=proj_map)
            diff_evapo_spatial.plot(transform=proj_data, cmap='RdBu_r', add_colorbar=True,
                                    subplot_kws={'projection': proj_map},
                                    cbar_kwargs={"label": "$\Delta$ ET [mm yr-1]", "shrink": 0.6})
            ax0.add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
            plt.title('')
            plt.axis('off')
            ax1 = plt.subplot(224, projection=proj_map)
            p = diff_swe_total_spatial.plot(transform=proj_data, cmap='RdBu_r', add_colorbar=True,
                                            subplot_kws={'projection': proj_map},
                                            cbar_kwargs={"label": "$\Delta$ SWE [mm yr-1]", "shrink": 0.6})
            ax1.add_feature(cf.BORDERS, linewidth=1, edgecolor='dimgray', alpha=1)
            plt.title('')
            plt.axis('off')
            plt.tight_layout()
            fig1.savefig(bf_plots_in / Path(comp_text+'_'+str(yr_in)+'_.png'),
                         facecolor='white', transparent=False)
            plt.close()
