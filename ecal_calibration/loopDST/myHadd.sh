#!/bin/bash
. /u/prozorov/wait_jobs.sh
filename=${1}
foldername=${2}
currentDir=$(pwd)


if [[ ${2} == *"_3200A"* ]]; then
 daystart=67
 dayfinish=89
fi

if [[ ${2} == *"_200A"* ]]; then
 daystart=90
 dayfinish=92
fi

if [[ ${2} == *"_2500A"* ]]; then
 daystart=88
 dayfinish=91
fi


for (( day=$daystart ; day<$dayfinish; day++ ))
do

    path2files="/lustre/hades/user/prozorov/mar19/leptonttree/gen5/${foldername}_0${day}"
    cd ${path2files}

nFiles=$( ls ${filename}_*.root | wc -l )
echo  "nFiles = ${nFiles}"

sOut="/lustre/hades/user/prozorov/mar19/leptonttree/gen5/${foldername}_${filename}_day0${day}leptons.root"
sAdd=""

for (( i=1 ; i<$nFiles ; i++ ))
do
sAdd="${sAdd} ${filename}_${i}.root"
done

hadd -f -k ${sOut} ${sAdd}
done

sOut="/lustre/hades/user/prozorov/mar19/leptonttree/gen5/all_days_leptonttree_gen5_${foldername}_${filename}.root"
filename="/lustre/hades/user/prozorov/mar19/leptonttree/gen5/${foldername}_${filename}_day0*leptons.root"
hadd -f -k ${sOut} ${filename}

cp  ${sOut} ${currentDir}/..

cd ${currentDir}
