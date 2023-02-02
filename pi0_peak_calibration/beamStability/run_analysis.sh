#!/bin/bash


currentDir=$(pwd)
path2files="/lustre/hebe/hades/user/prozorov/mar19/diphotons/pi0PeakCalibration/beamStability/"
filename="MixAll"
sOut="pi0PeakHists.root"
make
#. ~/wait_jobs.sh
. ./sendScript_SL.sh 0
. ~/wait_jobs.sh
cd ${path2files}
sAdd=""
rm ${sOut}
hadd -j -k -f ${sOut} ${filename}_*.root
cp  ${sOut} ${currentDir}
cd ${currentDir}
