#!/bin/bash

filename="AgAg165_emc_mom_cell_photon"
path2files="/lustre/hades/user/prozorov/prop17/ag1650ag_energy_correction"
currentDir=$(pwd)
outfile="photon_correction.root"
cd ${path2files}



rm ${filename}_29.root
rm ${filename}_31.root
hadd -f  ${outfile} ${filename}_*.root

cp ${outfile} ${currentDir}/../

cd ${currentDir}/../
root -l -n makeCorrection.C+

