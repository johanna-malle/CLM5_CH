#!/bin/sh

names_array=("CH-Aws" "CH-Cha" "CH-Dav" "CH-Fru" "CH-Lae" "CH-Oe2")

length=${#names_array[@]}
output=/home/malle/CLM5_install/scratch/archive/PTCLM_fluxnet_CRUplus05_new_realPFT/lnd/hist
CASE1=PTCLM_fluxnet_CRUplus05_new_realPFT

for ((i=0;i<=$length-1;i++)); do

new_name=${names_array[$i]}
old_name=CH-Aws

cp user_nl_clm_copy user_nl_clm1
cp user_datm.streams.txt.CLM1kmCRUJRA.Prec_copy.CLM_USRDAT user_datm.streams.txt.CLM1kmCRUJRA.Prec.CLM_USRDAT1
cp user_datm.streams.txt.CLM1kmCRUJRA.LP_copy.CLM_USRDAT user_datm.streams.txt.CLM1kmCRUJRA.LP.CLM_USRDAT1
cp user_datm.streams.txt.CLM1kmCRUJRA.Solr_copy.CLM_USRDAT user_datm.streams.txt.CLM1kmCRUJRA.Solr.CLM_USRDAT1
cp user_datm.streams.txt.CLM1kmCRUJRA.WZ_copy.CLM_USRDAT user_datm.streams.txt.CLM1kmCRUJRA.WZ.CLM_USRDAT1
cp user_datm.streams.txt.CLM1kmCRUJRA.TQ_copy.CLM_USRDAT user_datm.streams.txt.CLM1kmCRUJRA.TQ.CLM_USRDAT1
sed -e "s/$old_name/$new_name/g" user_nl_clm1 > user_nl_clm
sed -e "s/CH-Aws/$new_name/g" user_datm.streams.txt.CLM1kmCRUJRA.Prec.CLM_USRDAT1 > user_datm.streams.txt.CLM1kmCRUJRA.Prec.CLM_USRDAT
sed -e "s/CH-Aws/$new_name/g" user_datm.streams.txt.CLM1kmCRUJRA.TQ.CLM_USRDAT1 > user_datm.streams.txt.CLM1kmCRUJRA.TQ.CLM_USRDAT
sed -e "s/CH-Aws/$new_name/g" user_datm.streams.txt.CLM1kmCRUJRA.Solr.CLM_USRDAT1 > user_datm.streams.txt.CLM1kmCRUJRA.Solr.CLM_USRDAT
sed -e "s/CH-Aws/$new_name/g" user_datm.streams.txt.CLM1kmCRUJRA.WZ.CLM_USRDAT1 > user_datm.streams.txt.CLM1kmCRUJRA.WZ.CLM_USRDAT
sed -e "s/CH-Aws/$new_name/g" user_datm.streams.txt.CLM1kmCRUJRA.LP.CLM_USRDAT1 > user_datm.streams.txt.CLM1kmCRUJRA.LP.CLM_USRDAT
rm user_datm.streams.txt.CLM1kmCRUJRA.Prec.CLM_USRDAT1
rm user_datm.streams.txt.CLM1kmCRUJRA.TQ.CLM_USRDAT1
rm user_datm.streams.txt.CLM1kmCRUJRA.Solr.CLM_USRDAT1
rm user_datm.streams.txt.CLM1kmCRUJRA.WZ.CLM_USRDAT1
rm user_datm.streams.txt.CLM1kmCRUJRA.LP.CLM_USRDAT1
rm user_nl_clm1

./xmlchange ATM_DOMAIN_FILE=domain.lnd.${new_name}_navy.nc
./xmlchange LND_DOMAIN_FILE=domain.lnd.${new_name}_navy.nc

./case.build
./case.submit

for ren in $(seq 0 1 7)
do
mv $output/${CASE1}.clm2.h$ren.2003-01-01-00000.nc $output/${CASE1}_${new_name}.clm2.h$ren.2003-01-01-00000.nc
done

done