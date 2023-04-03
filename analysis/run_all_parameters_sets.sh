##!/usr/bin/env bash
currentDir=$(pwd)
path2files="/u/prozorov/mar19/diphotons/gen5a/systematic_uncertainty"
cd ${path2files}
rm nohup.out
array=($(ls -r parameters_files/.))
for i in "${array[@]}"
do
	  nohup bash run.sh ./parameters_files/${i} &
done
wait

 cd ana_results

root -l -n -q -b 'ana_results.cpp("../parameters_files/parameters_ANGLE_CUT1.conf ../parameters_files/parameters_REF0.conf ../parameters_files/parameters_RPC_TIME_CUT1.conf ../parameters_files/parameters_RPC_TIME_CUT2.conf ../parameters_files/parameters_N_MASS_BINS1.conf ../parameters_files/parameters_RPC_XY_CUT1.conf ../parameters_files/parameters_RPC_XY_CUT2.conf ../parameters_files/parameters_CLUSTER_SIZE_CUT2.conf ../parameters_files/parameters_ENERGY_CUT1.conf ../parameters_files/parameters_ENERGY_CUT2.conf")'
 root -l -n -q -b compareToChargedPions.cpp+
