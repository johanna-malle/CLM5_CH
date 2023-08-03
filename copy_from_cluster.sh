#!/bin/sh

bf=/media/malle/LaCie1/CLM5_SP_TILES/all_simulations_oldSurf
loc_hyp_archive=/home/malle/CESM2_install/scratch/archive

scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_03N_SP_CRUJRA_noLapse_OLD/lnd/hist/*.h* $bf/ch_1km_03N_SP_CRUJRA_noLapse/.
scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_03S_SP_CRUJRA_noLapse_OLD/lnd/hist/*.h* $bf/ch_1km_03S_SP_CRUJRA_noLapse/.

scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_02N_SP_CRUJRA_noLapse_OLD/lnd/hist/*.h* $bf/ch_1km_02N_SP_CRUJRA_noLapse/.
scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_02S_SP_CRUJRA_noLapse_OLD/lnd/hist/*.h* $bf/ch_1km_02S_SP_CRUJRA_noLapse/.

scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_04N_SP_CRUJRA_noLapse_OLD/lnd/hist/*.h* $bf/ch_1km_04N_SP_CRUJRA_noLapse/.
scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_04S_SP_CRUJRA_noLapse_OLD/lnd/hist/*.h* $bf/ch_1km_04S_SP_CRUJRA_noLapse/.

scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_01N_SP_OSHD/lnd/hist/*.h* $bf/ch_1km_01N_SP_OSHD/.
scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_01N_SP_CRUJRA/lnd/hist/*.h* $bf/ch_1km_01N_SP_CRUJRA/.

scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_01S_SP_OSHD/lnd/hist/*.h* $bf/ch_1km_01S_SP_OSHD/.
scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_01S_SP_CRUJRA/lnd/hist/*.h* $bf/ch_1km_01S_SP_CRUJRA/.

scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_05N_SP_OSHD/lnd/hist/*.h* $bf/ch_1km_05N_SP_OSHD/.
scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_05N_SP_CRUJRA/lnd/hist/*.h* $bf/ch_1km_05N_SP_CRUJRA/.

scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_05S_SP_OSHD/lnd/hist/*.h* $bf/ch_1km_05S_SP_OSHD/.
scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_05S_SP_CRUJRA/lnd/hist/*.h* $bf/ch_1km_05S_SP_CRUJRA/.

scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_04S_SP_OSHD/lnd/hist/*.h* $bf/ch_1km_04S_SP_OSHD/.

scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_1km_02N_SP_OSHD/lnd/hist/*.h* $bf/ch_1km_02N_SP_OSHD/.
scp malle@hyperion.wsl.ch:$loc_hyp_archive/ch_2km_02N_SP_CRUJRA/lnd/hist/*.h* $bf/ch_1km_02N_SP_CRUJRA/.