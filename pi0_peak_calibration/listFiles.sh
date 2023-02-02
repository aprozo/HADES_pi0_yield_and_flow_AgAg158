#!/bin/bash
#rm lists/photon_ttrees_gen5.list
currentDir=$(pwd)
outputdir=/lustre/hebe/hades/user/prozorov/mar19/photon_ttree_pi0_peak_corrected

for (( day_number=68; day_number<89;  day_number++ ))
do
  ls  ${outputdir}/day0${day_number}/photon_ttree_day0${day_number}_?.root >> lists/photon_ttrees_gen5_pi0_peak_corrected.list
  ls  ${outputdir}/day0${day_number}/photon_ttree_day0${day_number}_??.root >> lists/photon_ttrees_gen5_pi0_peak_corrected.list
  ls  ${outputdir}/day0${day_number}/photon_ttree_day0${day_number}_???.root >> lists/photon_ttrees_gen5_pi0_peak_corrected.list
done
