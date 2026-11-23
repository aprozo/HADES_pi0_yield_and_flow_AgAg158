#!/bin/bash


currentDir=$(pwd)
nInterations=15
path2files="/lustre/hades/user/prozorov/mar19/diphotons/pi0PeakCalibration"
filename="MixAll"
sOut="pi0PeakHists.root"
make
. ~/wait_jobs.sh

  for parN in 1.0
  do
  #root -l -n -b -q "extractCorrFactors.cpp+(0,${parN})"

    for (( iIter=0; iIter<$nInterations; iIter++ ))
      do
      # FILE=AllCorrFactorsIteration${iIter}parN${parN}.root
      # if  [ -f "$FILE" ]; then
      #   cp "$FILE" AllCorrFactors.root
      #   continue
      # fi

      echo par = ${parN} iteration = ${iIter}
      ( . ./sendScript_SL.sh ${iIter})
      echo par = ${parN} iteration = ${iIter}
      . ~/wait_jobs.sh
      cd ${path2files}
      sAdd=""
      rm ${sOut}
      hadd -j -k -f ${sOut} ${filename}_*.root
      cp  ${sOut} ${currentDir}
      rm  "${filename}_*.root"
      cd ${currentDir}
      echo par = ${parN} iteration = ${iIter}
      root -l -n -b -q "extractCorrFactors.cpp+(${iIter},${parN})"
      done

    #  cp AllCorrFactors.root AllCorrFactors${parN}.root
    #  rm AllCorrFactors.root

done
