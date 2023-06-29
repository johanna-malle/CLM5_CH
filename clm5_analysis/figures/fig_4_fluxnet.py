# -*- coding: utf-8 -*-
"""
Desc:
Created on 02.02.23 09:16
@author: malle
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import seaborn as sns
import statsmodels.formula.api as smf
from sklearn import preprocessing

def mae(y_true, predictions):
    y_true, predictions = np.array(y_true), np.array(predictions)
    return np.round(np.mean(np.abs(y_true - predictions)),2)

def rmse(y_true, predictions):
    y_true, predictions = np.array(y_true), np.array(predictions)
    return np.round(np.sqrt(((predictions - y_true) ** 2).mean()),2)

bf = Path('/home/lud11/malle/CLM5_CH/new_figures/comp_meas')

CH_Aws_abs = pd.read_csv(bf / 'CH-Aws_flux_comp.csv')
CH_Cha_abs = pd.read_csv(bf / 'CH-Cha_flux_comp.csv')
CH_Dav_abs = pd.read_csv(bf / 'CH-Dav_flux_comp.csv')
CH_Fru_abs = pd.read_csv(bf / 'CH-Fru_flux_comp.csv')
CH_Lae_abs = pd.read_csv(bf / 'CH-Lae_flux_comp.csv')
CH_Oe2_abs = pd.read_csv(bf / 'CH-Oe2_flux_comp.csv')

CH_Aws = pd.read_csv(bf / 'CH-Aws_flux_diff.csv')
CH_Cha = pd.read_csv(bf / 'CH-Cha_flux_diff.csv')
CH_Dav = pd.read_csv(bf / 'CH-Dav_flux_diff.csv')
CH_Fru = pd.read_csv(bf / 'CH-Fru_flux_diff.csv')
CH_Lae = pd.read_csv(bf / 'CH-Lae_flux_diff.csv')
CH_Oe2 = pd.read_csv(bf / 'CH-Oe2_flux_diff.csv')

df_Aws_lh = pd.melt(np.abs(CH_Aws[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']]),
                  var_name="model_run", value_name="lh")
df_Aws_lh['site']="CH-Aws"

df_Cha_lh = pd.melt(np.abs(CH_Cha[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']]),
                  var_name="model_run", value_name="lh")
df_Cha_lh['site']="CH-Cha"

df_Dav_lh = pd.melt(np.abs(CH_Dav[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']]),
                  var_name="model_run", value_name="lh")
df_Dav_lh['site']="CH-Dav"

df_Fru_lh = pd.melt(np.abs(CH_Fru[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']]),
                  var_name="model_run", value_name="lh")
df_Fru_lh['site']="CH-Fru"

df_Lae_lh = pd.melt(np.abs(CH_Lae[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']]),
                  var_name="model_run", value_name="lh")
df_Lae_lh['site']="CH-Lae"

df_Oe2_lh = pd.melt(np.abs(CH_Oe2[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']]),
                  var_name="model_run", value_name="lh")
df_Oe2_lh['site']="CH-Oe2"

frames = [df_Aws_lh, df_Cha_lh, df_Dav_lh,df_Fru_lh,df_Lae_lh,df_Oe2_lh]
combined_lh = pd.concat(frames)

# now sensible heat
df_Aws_sh = pd.melt(np.abs(CH_Aws[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']]),
                  var_name="model_run", value_name="sh")
df_Aws_sh['site']="CH-Aws"

df_Cha_sh = pd.melt(np.abs(CH_Cha[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']]),
                  var_name="model_run", value_name="sh")
df_Cha_sh['site']="CH-Cha"

df_Dav_sh = pd.melt(np.abs(CH_Dav[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']]),
                  var_name="model_run", value_name="sh")
df_Dav_sh['site']="CH-Dav"

df_Fru_sh = pd.melt(np.abs(CH_Fru[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']]),
                  var_name="model_run", value_name="sh")
df_Fru_sh['site']="CH-Fru"

df_Lae_sh = pd.melt(np.abs(CH_Lae[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']]),
                  var_name="model_run", value_name="sh")
df_Lae_sh['site']="CH-Lae"

df_Oe2_sh = pd.melt(np.abs(CH_Oe2[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']]),
                  var_name="model_run", value_name="sh")
df_Oe2_sh['site']="CH-Oe2"

frames = [df_Aws_sh, df_Cha_sh, df_Dav_sh,df_Fru_sh,df_Lae_sh,df_Oe2_sh]
combined_sh = pd.concat(frames)

# now gpp
df_Aws_gpp = pd.melt(np.abs(CH_Aws[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']]),
                  var_name="model_run", value_name="gpp")
df_Aws_gpp['site']="CH-Aws"

df_Cha_gpp = pd.melt(np.abs(CH_Cha[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']]),
                  var_name="model_run", value_name="gpp")
df_Cha_gpp['site']="CH-Cha"

df_Dav_gpp = pd.melt(np.abs(CH_Dav[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']]),
                  var_name="model_run", value_name="gpp")
df_Dav_gpp['site']="CH-Dav"

df_Fru_gpp = pd.melt(np.abs(CH_Fru[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']]),
                  var_name="model_run", value_name="gpp")
df_Fru_gpp['site']="CH-Fru"

df_Lae_gpp = pd.melt(np.abs(CH_Lae[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']]),
                  var_name="model_run", value_name="gpp")
df_Lae_gpp['site']="CH-Lae"

df_Oe2_gpp = pd.melt(np.abs(CH_Oe2[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']]),
                  var_name="model_run", value_name="gpp")
df_Oe2_gpp['site']="CH-Oe2"

frames = [df_Aws_gpp, df_Cha_gpp, df_Dav_gpp, df_Fru_gpp, df_Lae_gpp, df_Oe2_gpp]
combined_gpp = pd.concat(frames)

frames = [np.abs(CH_Aws[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']]),
         np.abs(CH_Cha[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']]),
         np.abs(CH_Dav[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']]),
         np.abs(CH_Fru[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']]),
         np.abs(CH_Lae[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']]),
         np.abs(CH_Oe2[['cru_05_lh','cru_1km_lh','cru_plus_05_lh','cru_plus_1km_lh','oshd_05_lh','oshd_1km_lh']])]
dh_lh = pd.concat(frames)

frames = [np.abs(CH_Aws[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']]),
         np.abs(CH_Cha[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']]),
         np.abs(CH_Dav[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']]),
         np.abs(CH_Fru[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']]),
         np.abs(CH_Lae[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']]),
         np.abs(CH_Oe2[['cru_05_sh','cru_1km_sh','cru_plus_05_sh','cru_plus_1km_sh','oshd_05_sh','oshd_1km_sh']])]
dh_sh = pd.concat(frames)

frames = [np.abs(CH_Aws[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']]),
         np.abs(CH_Cha[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']]),
         np.abs(CH_Dav[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']]),
         np.abs(CH_Fru[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']]),
         np.abs(CH_Lae[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']]),
         np.abs(CH_Oe2[['cru_05_gpp','cru_1km_gpp','cru_plus_05_gpp','cru_plus_1km_gpp','oshd_05_gpp','oshd_1km_gpp']])]
dh_gpp = pd.concat(frames)

frames = [np.abs(CH_Aws[['oshd_1km_h2osoi','oshd_05_h2osoi','cru_1km_h2osoi','cru_05_h2osoi','cru_plus_1km_h2osoi','cru_plus_05_h2osoi']]),
         np.abs(CH_Cha[['oshd_1km_h2osoi', 'oshd_05_h2osoi','cru_1km_h2osoi','cru_05_h2osoi','cru_plus_1km_h2osoi','cru_plus_05_h2osoi']]),
         np.abs(CH_Dav[['oshd_1km_h2osoi','oshd_05_h2osoi','cru_1km_h2osoi','cru_05_h2osoi','cru_plus_1km_h2osoi','cru_plus_05_h2osoi']]),
         np.abs(CH_Fru[['oshd_1km_h2osoi','oshd_05_h2osoi','cru_1km_h2osoi','cru_05_h2osoi','cru_plus_1km_h2osoi','cru_plus_05_h2osoi']]),
         np.abs(CH_Lae[['oshd_1km_h2osoi','oshd_05_h2osoi','cru_1km_h2osoi','cru_05_h2osoi','cru_plus_1km_h2osoi','cru_plus_05_h2osoi']]),
         np.abs(CH_Oe2[['oshd_1km_h2osoi','oshd_05_h2osoi','cru_1km_h2osoi','cru_05_h2osoi','cru_plus_1km_h2osoi','cru_plus_05_h2osoi']])]
dh_swc = pd.concat(frames)

my_pal_lh = {"oshd_1km_lh":(27/255,158/255,119/255), "oshd_05_lh": (27/255,158/255,119/255),
             "cru_plus_1km_lh":(117/255,112/255,179/255), "cru_plus_05_lh":(117/255,112/255,179/255),
             "cru_1km_lh":(217/255,95/255,2/255), "cru_05_lh":(217/255,95/255,2/255)}

my_pal_sh = {"oshd_1km_sh": (27/255,158/255,119/255), "oshd_05_sh": (27/255,158/255,119/255),
             "cru_plus_1km_sh":(117/255,112/255,179/255), "cru_plus_05_sh":(117/255,112/255,179/255),
             "cru_1km_sh":(217/255,95/255,2/255), "cru_05_sh":(217/255,95/255,2/255)}

my_pal_gpp = {"oshd_1km_gpp": (27/255,158/255,119/255), "oshd_05_gpp": (27/255,158/255,119/255),
              "cru_plus_1km_gpp":(117/255,112/255,179/255), "cru_plus_05_gpp":(117/255,112/255,179/255),
              "cru_1km_gpp":(217/255,95/255,2/255), "cru_05_gpp":(217/255,95/255,2/255)}


# fit linear mixed effects mdoel:
# r-code: lm1 = lmer(TSS ~ group + (1|id) + (1|model),data=dfx,REML=TRUE)

gpp_act=combined_gpp[~np.isnan(combined_gpp.gpp)].reset_index(drop=True)
sh_act=combined_sh[~np.isnan(combined_sh.sh)].reset_index(drop=True)
lh_act=combined_lh[~np.isnan(combined_lh.lh)].reset_index(drop=True)

rng = (0, 1)  # specify the range to which you want to scale
scaler = preprocessing.MinMaxScaler(feature_range=(rng[0], rng[1]))
normed_gpp = scaler.fit_transform(np.array(gpp_act.gpp).reshape(-1, 1))
normed_lh = scaler.fit_transform(np.array(lh_act.lh).reshape(-1, 1))
normed_sh = scaler.fit_transform(np.array(sh_act.sh).reshape(-1, 1))

gpp_act_norm = gpp_act.copy()
gpp_act_norm['gpp']=normed_gpp
sh_act_norm = sh_act.copy()
sh_act_norm['sh']=normed_sh
lh_act_norm = lh_act.copy()
lh_act_norm['lh']=normed_lh
#norm_lst = [round(i[0],2) for i in normed]  # the output is an array of arrays, so tidy the dimensions

lme_out_gpp = smf.mixedlm('gpp ~ model_run', gpp_act, groups=gpp_act['site'])
lme_fit_gpp = lme_out_gpp.fit()
lme_out_gpp_norm = smf.mixedlm('gpp ~ model_run', gpp_act_norm, groups=gpp_act_norm['site'])
lme_fit_gpp_norm = lme_out_gpp_norm.fit()

lme_out_sh = smf.mixedlm('sh ~ model_run', sh_act, groups=sh_act['site'])
lme_fit_sh = lme_out_sh.fit()
lme_out_sh_norm = smf.mixedlm('sh ~ model_run', sh_act_norm, groups=sh_act_norm['site'])
lme_fit_sh_norm = lme_out_sh_norm.fit()

lme_out_lh = smf.mixedlm('lh ~ model_run', lh_act, groups=lh_act['site'])
lme_fit_lh = lme_out_lh.fit()
lme_out_lh_norm = smf.mixedlm('lh ~ model_run', lh_act_norm, groups=lh_act_norm['site'])
lme_fit_lh_norm = lme_out_lh_norm.fit()

flierprops = dict(marker='o', markersize=4.5, markeredgecolor='gray', markerfacecolor='silver', alpha=0.45)
color_all = [(217/255,95/255,2/255),(117/255,112/255,179/255),(117/255,112/255,179/255),
             (27/255,158/255,119/255),(27/255,158/255,119/255)]
# make comparative boxplots
# creating grid for subplots
fig = plt.figure(layout="constrained", figsize=(14, 6))
ax1 = plt.subplot2grid(shape=(5, 1), loc=(0, 0), rowspan=4)
ax=sns.boxplot(x="site", y="lh", hue="model_run", palette=my_pal_lh, data=combined_lh, flierprops=flierprops, showmeans=True,
               meanprops={"marker":"o", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"4"},
               medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), linewidth=1.6, saturation=0.8, showfliers = False)
sns.despine(offset=10, trim=True)
plt.ylabel('($\Delta$ LH)$_{abs}$ [W m-2]',fontsize=17)
plt.xlabel('')

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(15)
for patch in ax1.patches[1::2]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .5))
for patch in ax1.patches[1:10:4]: #1
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 1))

le1=ax1.legend(loc='upper left', title="", ncols=3, bbox_to_anchor=(0.28, 0.98),fontsize=13)
le1.texts[0].set_text('Clim$_{CRU 0.5^{\circ}}$')
le1.texts[1].set_text('Clim$_{CRU 1km}$')
le1.texts[2].set_text('Clim$_{CRU* 0.5^{\circ}}$')
le1.texts[3].set_text('Clim$_{CRU* 1km}$')
le1.texts[4].set_text('Clim$_{OSHD 0.5^{\circ}}$')
le1.texts[5].set_text('Clim$_{OSHD 1km}$')
for lh in le1.legend_handles[1::2]:
    lh.set_alpha(0.5)
ti1 = ax1.set_title('(a)', loc='left', x=-0.09, y=1.03, pad=-24, fontsize=18, fontweight='bold') #,

ax1a = plt.subplot2grid(shape=(5, 1), loc=(4, 0), rowspan=1)
a1=lme_fit_lh.summary().tables[1]
st_1=a1['[0.025'][1:-1].values.astype('float')
st_2=a1['0.975]'][1:-1].values.astype('float')
a=np.asarray([[[st_1[0],1],[st_2[0],1]],[[st_1[1],2],[st_2[1],2]],[[st_1[2],3],[st_2[2],3]],[[st_1[3],4],[st_2[3],4]],
              [[st_1[4],5],[st_2[4],5]]])
X=a[:,:,0].T
Y=a[:,:,1].T

for i in range(1, 6):
    if i in [1,3,5]:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=8, alpha=0.5)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1], alpha=0.5, linewidth=2)
    else:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=8)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1], linewidth=2)
plt.axvline(x=0, linewidth=2.5, linestyle='--', color='r')
ax1a.set_yticks([1,2,3,4,5])
ax1a.set_yticklabels(['Clim$_{CRU 1km}$','Clim$_{CRU* 0.5^{\circ}}$','Clim$_{CRU* 1km}$','Clim$_{OSHD 0.5^{\circ}}$',
                      'Clim$_{OSHD 1km}$'],fontsize=13.5)
ax1a.set_xlabel('Coefficient estimate',fontsize=16)
ax1a.set_ylim([0,6])
for label in ax1a.get_xticklabels():
	label.set_fontsize(15)

plt.show()
fig.savefig(bf / 'figure_4_a.png', facecolor='white', transparent=False)
fig.savefig(bf / 'figure_4_a.pdf', facecolor='white', transparent=False)


# now sensible heat
fig = plt.figure(layout="constrained", figsize=(14, 6))
ax2 = plt.subplot2grid(shape=(5, 1), loc=(0, 0), rowspan=4)
ax=sns.boxplot(x="site", y="sh", hue="model_run", palette=my_pal_sh, data=combined_sh, flierprops=flierprops, showmeans=True,
               meanprops={"marker":"o", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"4"},
               medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'),
               linewidth=1.6, saturation=0.8, showfliers = False)
ax2.get_legend().remove()
sns.despine(offset=10, trim=True)
plt.ylabel('($\Delta$ SH)$_{abs}$ [W m-2]',fontsize=17)
plt.xlabel('')

for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
	label.set_fontsize(15)

for patch in ax2.patches[1::2]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .5))

for patch in ax2.patches[1:10:4]: #1
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 1))

ti1 = ax2.set_title('(b)', loc='left', x=-0.09, y=1.1, pad=-24, fontsize=18, fontweight='bold')

ax2a = plt.subplot2grid(shape=(5, 1), loc=(4, 0), rowspan=1)
a1=lme_fit_sh.summary().tables[1]
st_1=a1['[0.025'][1:-1].values.astype('float')
st_2=a1['0.975]'][1:-1].values.astype('float')
a=np.asarray([[[st_1[0],1],[st_2[0],1]],[[st_1[1],2],[st_2[1],2]],[[st_1[2],3],[st_2[2],3]],[[st_1[3],4],[st_2[3],4]],
              [[st_1[4],5],[st_2[4],5]]])
X=a[:,:,0].T
Y=a[:,:,1].T
for i in range(1, 6):
    if i in [1,3,5]:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=8, alpha=0.5)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1], alpha=0.5, linewidth=2)
    else:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=8)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1], linewidth=2)

plt.axvline(x=0, linewidth=2.5, linestyle='--', color='r')
ax2a.set_yticks([1,2,3,4,5])
ax2a.set_yticklabels(['Clim$_{CRU 1km}$','Clim$_{CRU* 0.5^{\circ}}$','Clim$_{CRU* 1km}$','Clim$_{OSHD 0.5^{\circ}}$',
                      'Clim$_{OSHD 1km}$'],fontsize=13.5)
ax2a.set_xlabel('Coefficient estimate',fontsize=16)
ax2a.set_ylim([0, 6])
for label in ax2a.get_xticklabels():
	label.set_fontsize(15)

plt.show()
fig.savefig(bf / 'figure_4_b.png', facecolor='white', transparent=False)
fig.savefig(bf / 'figure_4_b.pdf', facecolor='white', transparent=False)

# now gpp
fig = plt.figure(layout="constrained", figsize=(14, 6))
ax3 = plt.subplot2grid(shape=(5, 1), loc=(0, 0), rowspan=4)
ax=sns.boxplot(x="site", y="gpp", hue="model_run", palette=my_pal_gpp, data=combined_gpp, flierprops=flierprops, showmeans=True,
               meanprops={"marker":"o", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"4"},
               medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), linewidth=1.6, saturation=0.8, showfliers = False)
ax3.get_legend().remove()
sns.despine(offset=10, trim=True)
plt.ylabel('($\Delta$ GPP)$_{abs}$ [umol m-2 s-1]', fontsize=17)
plt.xlabel('')

for label in (ax3.get_xticklabels() + ax3.get_yticklabels()):
	label.set_fontsize(15)

for patch in ax3.patches[1::2]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .5))

for patch in ax3.patches[1:10:4]: #1
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 1))

ti1 = ax3.set_title('(c)', loc='left', x=-0.09, y=1.03, pad=-24, fontsize=18,fontweight='bold') #,

ax3a = plt.subplot2grid(shape=(5, 1), loc=(4, 0), rowspan=1)
a1=lme_fit_gpp.summary().tables[1]
st_1=a1['[0.025'][1:-1].values.astype('float')
st_2=a1['0.975]'][1:-1].values.astype('float')
a=np.asarray([[[st_1[0],1],[st_2[0],1]],[[st_1[1],2],[st_2[1],2]],[[st_1[2],3],[st_2[2],3]],[[st_1[3],4],[st_2[3],4]],
              [[st_1[4],5],[st_2[4],5]]])
X=a[:,:,0].T
Y=a[:,:,1].T

for i in range(1, 6):
    if i in [1,3,5]:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=8, alpha=0.5)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1], alpha=0.5, linewidth=2)
    else:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=8)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1], linewidth=2)

plt.axvline(x=0, linewidth=2.5, linestyle='--', color='r')
ax3a.set_yticks([1,2,3,4,5])
ax3a.set_yticklabels(['Clim$_{CRU 1km}$','Clim$_{CRU* 0.5^{\circ}}$','Clim$_{CRU* 1km}$','Clim$_{OSHD 0.5^{\circ}}$',
                      'Clim$_{OSHD 1km}$'], fontsize=13.5)
ax3a.set_xlabel('Coefficient estimate', fontsize=16)
ax3a.set_ylim([0,6])
for label in ax3a.get_xticklabels():
	label.set_fontsize(15)
plt.show()
fig.savefig(bf / 'figure_4_c.png', facecolor='white', transparent=False)
fig.savefig(bf / 'figure_4_c.pdf', facecolor='white', transparent=False)


#### now same plot with normed data
fig = plt.figure(layout="constrained", figsize=(14, 5))
ax1 = plt.subplot2grid(shape=(5, 1), loc=(0, 0), rowspan=4)
ax=sns.boxplot(x="site", y="lh", hue="model_run", palette=my_pal_lh, data=combined_lh, flierprops=flierprops, showmeans=True,
               meanprops={"marker":"o", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"4"},
               medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), linewidth=1.6, saturation=0.8, showfliers = False)
sns.despine(offset=10, trim=True)
plt.ylabel('($\Delta$ LH)$_{abs}$ [W m-2]',fontsize=13)
plt.xlabel('')

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(12)
for patch in ax1.patches[1::2]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .5))
for patch in ax1.patches[1:10:4]: #1
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 1))

le1=ax1.legend(loc='upper left', title="", ncols=3, bbox_to_anchor=(0.28, 0.98),fontsize=14)
le1.texts[0].set_text('Clim$_{CRU 0.5^{\circ}}$')
le1.texts[1].set_text('Clim$_{CRU 1km}$')
le1.texts[2].set_text('Clim$_{CRU* 0.5^{\circ}}$')
le1.texts[3].set_text('Clim$_{CRU* 1km}$')
le1.texts[4].set_text('Clim$_{OSHD 0.5^{\circ}}$')
le1.texts[5].set_text('Clim$_{OSHD 1km}$')
for lh in le1.legend_handles[1::2]:
    lh.set_alpha(0.5)
ti1 = ax1.set_title('(a)', loc='left', x=-0.05, y=1.03, pad=-24, fontsize=15, fontweight='bold') #,

ax1a = plt.subplot2grid(shape=(5, 1), loc=(4, 0), rowspan=1)
a1=lme_fit_lh_norm.summary().tables[1]
st_1=a1['[0.025'][1:-1].values.astype('float')
st_2=a1['0.975]'][1:-1].values.astype('float')
a=np.asarray([[[st_1[0],1],[st_2[0],1]],[[st_1[1],2],[st_2[1],2]],[[st_1[2],3],[st_2[2],3]],[[st_1[3],4],[st_2[3],4]],
              [[st_1[4],5],[st_2[4],5]]])
X=a[:,:,0].T
Y=a[:,:,1].T

for i in range(1, 6):
    if i in [1,3,5]:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=6, alpha=0.5)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1], alpha=0.5, linewidth=1.4)
    else:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=6)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1], linewidth=1.4)
plt.axvline(x=0, linewidth=2, linestyle='--', color='r')
ax1a.set_yticks([1,2,3,4,5])
ax1a.set_yticklabels(['Clim$_{CRU 1km}$','Clim$_{CRU* 0.5^{\circ}}$','Clim$_{CRU* 1km}$','Clim$_{OSHD 0.5^{\circ}}$',
                      'Clim$_{OSHD 1km}$'],fontsize=10.5)
ax1a.set_xlabel('Coefficient estimate',fontsize=13)
ax1a.set_ylim([0,6])
for label in ax1a.get_xticklabels():
	label.set_fontsize(12)

plt.show()
fig.savefig(bf / 'figure_4_a_norm.png', facecolor='white', transparent=False)
fig.savefig(bf / 'figure_4_a_norm.pdf', facecolor='white', transparent=False)


# now sensible heat
fig = plt.figure(layout="constrained", figsize=(14, 5))
ax2 = plt.subplot2grid(shape=(5, 1), loc=(0, 0), rowspan=4)
ax=sns.boxplot(x="site", y="sh", hue="model_run", palette=my_pal_sh, data=combined_sh, flierprops=flierprops, showmeans=True,
               meanprops={"marker":"o", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"4"},
               medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'),
               linewidth=1.6, saturation=0.8, showfliers = False)
ax2.get_legend().remove()
sns.despine(offset=10, trim=True)
plt.ylabel('($\Delta$ SH)$_{abs}$ [W m-2]',fontsize=13)
plt.xlabel('')

for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
	label.set_fontsize(12)

for patch in ax2.patches[1::2]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .5))

for patch in ax2.patches[1:10:4]: #1
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 1))

ti1 = ax2.set_title('(b)', loc='left', x=-0.08, y=1.1, pad=-24, fontsize=15, fontweight='bold')

ax2a = plt.subplot2grid(shape=(5, 1), loc=(4, 0), rowspan=1)
a1=lme_fit_sh_norm.summary().tables[1]
st_1=a1['[0.025'][1:-1].values.astype('float')
st_2=a1['0.975]'][1:-1].values.astype('float')
a=np.asarray([[[st_1[0],1],[st_2[0],1]],[[st_1[1],2],[st_2[1],2]],[[st_1[2],3],[st_2[2],3]],[[st_1[3],4],[st_2[3],4]],
              [[st_1[4],5],[st_2[4],5]]])
X=a[:,:,0].T
Y=a[:,:,1].T
for i in range(1, 6):
    if i in [1,3,5]:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=4, alpha=0.5)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1], alpha=0.5)
    else:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=4)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1])

plt.axvline(x=0, linewidth=2, linestyle='--', color='r')
ax2a.set_yticks([1,2,3,4,5])
ax2a.set_yticklabels(['Clim$_{CRU 1km}$','Clim$_{CRU* 0.5^{\circ}}$','Clim$_{CRU* 1km}$','Clim$_{OSHD 0.5^{\circ}}$',
                      'Clim$_{OSHD 1km}$'],fontsize=10.5)
ax2a.set_xlabel('Coefficient estimate',fontsize=13)
ax2a.set_ylim([0, 6])
for label in ax2a.get_xticklabels():
	label.set_fontsize(12)

plt.show()
fig.savefig(bf / 'figure_4_b_norm.png', facecolor='white', transparent=False)

# now gpp
fig = plt.figure(layout="constrained", figsize=(14, 5))
ax3 = plt.subplot2grid(shape=(5, 1), loc=(0, 0), rowspan=4)
ax=sns.boxplot(x="site", y="gpp", hue="model_run", palette=my_pal_gpp, data=combined_gpp, flierprops=flierprops, showmeans=True,
               meanprops={"marker":"o", "markerfacecolor":"white", "markeredgecolor":"black", "markersize":"4"},
               medianprops=dict(color="grey", alpha=0.85, linewidth=1.9, linestyle='-'), linewidth=1.6, saturation=0.8, showfliers = False)
ax3.get_legend().remove()
sns.despine(offset=10, trim=True)
plt.ylabel('($\Delta$ GPP)$_{abs}$ [umol m-2 s-1]', fontsize=13)
plt.xlabel('')

for label in (ax3.get_xticklabels() + ax3.get_yticklabels()):
	label.set_fontsize(12)

for patch in ax3.patches[1::2]:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .5))

for patch in ax3.patches[1:10:4]: #1
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 1))

ti1 = ax3.set_title('(c)', loc='left', x=-0.08, y=1.03, pad=-24, fontsize=15,fontweight='bold') #,

ax3a = plt.subplot2grid(shape=(5, 1), loc=(4, 0), rowspan=1)
a1=lme_fit_gpp_norm.summary().tables[1]
st_1=a1['[0.025'][1:-1].values.astype('float')
st_2=a1['0.975]'][1:-1].values.astype('float')
a=np.asarray([[[st_1[0],1],[st_2[0],1]],[[st_1[1],2],[st_2[1],2]],[[st_1[2],3],[st_2[2],3]],[[st_1[3],4],[st_2[3],4]],
              [[st_1[4],5],[st_2[4],5]]])
X=a[:,:,0].T
Y=a[:,:,1].T

for i in range(1, 6):
    if i in [1,3,5]:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=4, alpha=0.5)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1], alpha=0.5)
    else:
        plt.plot(np.array(a1['Coef.'].values[i]).astype('float'),i,'o','filled',color=color_all[i-1],markersize=4)
        plt.plot(X[:,i-1], Y[:,i-1], color=color_all[i-1])

plt.axvline(x=0, linewidth=2, linestyle='--', color='r')
ax3a.set_yticks([1,2,3,4,5])
ax3a.set_yticklabels(['Clim$_{CRU 1km}$','Clim$_{CRU* 0.5^{\circ}}$','Clim$_{CRU* 1km}$','Clim$_{OSHD 0.5^{\circ}}$',
                      'Clim$_{OSHD 1km}$'], fontsize=10.5)
ax3a.set_xlabel('Coefficient estimate', fontsize=13)
ax3a.set_ylim([0,6])
for label in ax3a.get_xticklabels():
	label.set_fontsize(12)
plt.show()
fig.savefig(bf / 'figure_4_c_norm.png', facecolor='white', transparent=False)


