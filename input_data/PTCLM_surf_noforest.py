# -*- coding: utf-8 -*-
"""
Desc:
Created on 25.01.23 15:23
@author: malle
"""

import numpy as np
from pathlib import Path
import xarray as xr
import sys
import platform


if platform.system() == 'Windows':
    base_dir = Path('L:\malle\CLM5_CH')
    str_cut = '\\'
else:
    base_dir = Path('/home/lud11/malle/CLM5_CH')
    str_cut = '/'

ptclm5_old_orig_in = base_dir / 'analysis_dec_2021' / 'PTCLM_old_surfdata_orig'
ptclm5_old_in = base_dir / 'analysis_dec_2021' / 'PTCLM_old_surfdata'
ptclm5_new_in = base_dir / 'analysis_dec_2021' / 'PTCLM_new_surfdata'

ptclm5_old_orig_out = base_dir / 'analysis_dec_2021' / 'PTCLM_old_surfdata_orig_nofor'
ptclm5_old_out = base_dir / 'analysis_dec_2021' / 'PTCLM_old_surfdata_nofor'
ptclm5_new_out = base_dir / 'analysis_dec_2021' / 'PTCLM_new_surfdata_nofor'

files_in = [ptclm5_old_orig_in, ptclm5_old_in, ptclm5_new_in]
files_out = [ptclm5_old_orig_out, ptclm5_old_out, ptclm5_new_out]

for num_in, path_in in enumerate(files_in):
    path_out = files_out[num_in]

    files = path_in.glob('**/*')
    files1 = [x for x in files if x.is_file()]

    for file_in in files1:
        name_out = str(file_in).split(str_cut)[-1]
        file_out = path_out / Path(name_out)

        surf_in = xr.open_dataset(file_in)

        surf_in.LANDFRAC_PFT.values[:] = 1
        surf_in.PCT_GLACIER.values[:] = 0
        surf_in.PCT_LAKE.values[:] = 0
        surf_in.PCT_WETLAND.values[:] = 0
        surf_in.PCT_URBAN.values[:] = 0
        surf_in.PCT_NATVEG.values[:] = 0
        surf_in.PCT_CROP.values[:] = 0

        surf_in.PCT_NATVEG.values[:] = 100
        surf_in.PCT_NAT_PFT.values[:] = 0
        surf_in.PCT_NAT_PFT.values[0][:] = 100  # 100% bare ground

        surf_in.to_netcdf(file_out)
