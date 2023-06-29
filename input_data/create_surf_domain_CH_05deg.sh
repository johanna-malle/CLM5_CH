#!/bin/bash

#SBATCH -J regrid_all       # job name
#SBATCH -N 1
## SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00 #batch job time limit
#SBATCH -p 24hour
#SBATCH -o regrid_all.%J.out   # output filename
#SBATCH -e regrid_all.%J.err   # error filename

#----------------------------------------------------------------------
# Set parameters
#----------------------------------------------------------------------

names_array=("CH_1deg_v2")
lat_array=("46.845283509")
long_array=("8.172828382")
length=${#names_array[@]}

for ((i=0;i<=$length-1;i++)); do

export PTNAME=${names_array[$i]}
export p_LON=${long_array[$i]}
export p_LAT=${lat_array[$i]}
export NX="10"
export NY="6"
export CDATE="210804"
export ESMFBIN_PATH="/home/malle/ESMF/install_ESMF_8_1/bin"
export CSMDATA="/media/malle/LaCie1/Oswald/input_data"
export REGRID_PROC="1"

cd /home/malle/CLM5_install/my_cesm_sandbox/components/clm/tools/mkmapdata/
./mknoocnmap.pl -n $PTNAME -nx $NX -ny $NY -p $p_LAT,$p_LON -dy 3 -dx 5

export MAPFILE=/home/malle/CLM5_install/my_cesm_sandbox/components/clm/tools/mkmapgrids/map_${PTNAME}_noocean_to_${PTNAME}_nomask_aave_da_${CDATE}.nc
export GRIDFILE=/home/malle/CLM5_install/my_cesm_sandbox/components/clm/tools/mkmapgrids/SCRIPgrid_${PTNAME}_nomask_c${CDATE}.nc
export RES=1km-merge-10min

if [ -z "$RES" ]; then
   echo "Run for all valid resolutions"
   resols=`../../bld/queryDefaultNamelist.pl -res list -silent`
else
   resols="$RES"
fi
echo "Create mapping files for this list of resolutions: $resols"

#----------------------------------------------------------------------

for res in $resols; do
   echo "Create mapping files for: $res"
#----------------------------------------------------------------------
   cmdargs="-r $PTNAME"
   res_type="regional"
   cmdargs="$cmdargs -t $res_type -f $GRIDFILE -v"

   echo "$res_type"
       regrid_num_proc=$SLURM_NTASKS
       if [ ! -z "$SLURM_JOB_ACCOUNT" ]; then
           echo "batch"
	   cmdargs="$cmdargs -b"
       fi
   echo "args: $cmdargs"
   echo "time env REGRID_PROC=$regrid_num_proc ./mkmapdata.sh $cmdargs\n"
   time env REGRID_PROC=$regrid_num_proc ./mkmapdata.sh $cmdargs
done


cd /home/malle/CLM5_install/my_cesm_sandbox/cime/tools/mapping/gen_domain_files/
ulimit -s unlimited

export gridlnd=$PTNAME
export gridocn="navy"

./gen_domain -m $MAPFILE -o $gridocn -l $gridlnd

cd /home/malle/CLM5_install/my_cesm_sandbox/components/clm/tools/mksurfdata_map/
export HDF5_USE_FILE_LOCKING=FALSE
./mksurfdata.pl -r usrspec -usr_gname $PTNAME -usr_gdate $CDATE #-hirespft -glc -allownofile


done
