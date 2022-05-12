#!/bin/bash

currentDir=$(pwd)
SL_mon.pl>out.txt
grepout=$(grep -o 'prozorov' 'out.txt')

while [ ! -z "$grepout" ]
do
   echo "sleep another minute"
   echo $grepout
   sleep 30;
   SL_mon.pl>out.txt
   grepout=$(grep -o 'prozorov' 'out.txt')
done

filename="exp_lepton_ttree_gen3_29Jul"
path2files="/lustre/hebe/hades/user/prozorov/mar19/exp_lepton_ttree"


for (( day=67; day<=88; day++ ))
do
    cd ${path2files}/day0${day}
    nFiles=$( ls ${filename}_*.root | wc -l )
    echo  "nFiles = ${nFiles}"
    sAdd=""
    sOut="date$(date +"%d%b")_day${day}.root"
    for (( j=1; j<=$nFiles; j++ ))
        do
	    sAdd="${sAdd} ${filename}_${j}.root"
        done
    hadd -f -k ${sOut} ${sAdd}
    cp  ${sOut} ../
done

cd ${path2files}
output_ttree="exp_leptons_ttree_all_days.root"
input_ttrees="date$(date +"%d%b")_day*.root"
hadd -f -k ${output_ttree} ${input_ttrees}
cp ${output_ttree} ${currentDir}/../makeplots/
cd ${currentDir}
