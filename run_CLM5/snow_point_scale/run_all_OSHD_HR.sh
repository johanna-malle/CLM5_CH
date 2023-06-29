#!/bin/sh

names_array=("ZER3" "GUG2" "TRU2" "JUL2" "CMA2" "FRA2" "OBM2" "ILI2" "CBS" "ALT" "FRI" "BSG" "BLZ" "MAG" "1LC" "4MS" "5DF" "5DO" "5LQ" "5WJ" "6SB" "7BR" "7SC" "7SD" "7SU" "7ZN" "7ZU" "ABG" "AIR" "APT" "ARO" "CEO" "DAV2" "DAV3" "DAV4" "DAV5" "DEH" "DIS" "GOR2" "LAG2" "LAU2" "MAE2" "MAS" "OFE2" "SLF2" "SON" "VAL2" "VLS2" "WHA" "YBR2" "ZNZ2")

length=${#names_array[@]}
run=/home/malle/CLM5_install/scratch/PTCLM_all_OSHD_newSurf/run
output=/home/malle/CLM5_install/scratch/archive/PTCLM_all_OSHD_newSurf/lnd/hist
rest=/home/malle/CLM5_install/scratch/archive/PTCLM_all_crujra_lapse_spinup/rest
CASE1=PTCLM_all_OSHD_newSurf

for ((i=0;i<=$length-1;i++)); do

new_name=${names_array[$i]}
old_name=1LC

cp user_nl_clm_copy user_nl_clm1
cp user_datm.streams.txt.CLM1PT_copy.CLM_USRDAT user_datm.streams.txt.CLM1PT.CLM_USRDAT1

cp $rest/${new_name}/* $run/.

sed -e "s/$old_name/$new_name/g" user_nl_clm1 > user_nl_clm
sed -e "s/1LC/$new_name/g" user_datm.streams.txt.CLM1PT.CLM_USRDAT1 > user_datm.streams.txt.CLM1PT.CLM_USRDAT

rm user_datm.streams.txt.CLM1PT.CLM_USRDAT1
rm user_nl_clm1

./xmlchange ATM_DOMAIN_FILE=domain.lnd.${new_name}_navy.nc
./xmlchange LND_DOMAIN_FILE=domain.lnd.${new_name}_navy.nc

./case.build
./case.submit

for ren in $(seq 0 1 7)
do
mv $output/${CASE1}.clm2.h$ren.2014-01-01-00000.nc $output/${CASE1}_${new_name}.clm2.h$ren.2014-01-01-00000.nc
done

done
