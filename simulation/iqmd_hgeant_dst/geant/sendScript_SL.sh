#!/bin/bash

# Submission script for GridEngine (GE). Each job will 
# be executed via the jobScript.sh
# This jobScript supports up to 7 parameters. Edit 
# the user specific part of the script according to 
# your program.
#
# Input to the script is a filelist with 1 file per line.
# For each file a job is started. With the parameter 
# nFilesPerJob a comma separated filelist will be 
# generated and handed to the job script. This feature
# is usefull when running many small jobs. Each
# job has its own logfile. All needed directories for the
# logfiles will be created if non existing.
#
# IMPORTANT: the hera/prometheus cluster jobs will only
# see the /hera file system. All needed scripts, programs
# and parameters have to be located on /hera or the job
# will crash. This script syncs your working dir to the submission
# dir on /hera . Make sure your scripts use the submission dir!
# Software should be taken from /cvmfs/hades.gsi.de/install/
#
# job log files will be named like inputfiles. If nFilesPerJob > 1
# the log files will contain the partnumber.
#
######################################################################
#   CONFIGURATION

user=$(whoami)
currentDir=$(pwd | xargs -i basename {})
currentDir=../$currentDir

# Au+Au
impact=iqmd
#impact=bmax16

nevents=1000
nFiles=10000
         
#part=phi # minbias bmax10
part=iqmd # minbias bmax10
offset=0

submmissionbase=/lustre/hades/user/${user}/submit/flowtest
submissiondir=${submmissionbase}/geant
 nFilesPerJob=1                            # number of files to be analyzed by 1 job (default==1)
    jobscript=${submissiondir}/jobScript_SL.sh  # exec script (full path, call without dot, set it executable!)
    outputdir=/lustre/hades/user/${user}/flowtest/geant/     # outputdir for files AND logFiles
pathoutputlog=${outputdir}/log                  # protocol from batch farm for each file
     filename=testrun                           # filename of log file if nFilesPerJob > 1 (partnumber will be appended)
     replace=${submissiondir}/replaceIniDat.pl
     template=${submissiondir}/au123au.dat     # standard au123au.dat (random vertex,eventplane rot)  pluto:  au123au_pluto.dat (vertex from file, not eventplane rot)
      input1=/lustre/hades/user/harabasz/flowtest/evt/iqmd_evt      # bmax10
      output=/lustre/hades/user/harabasz/flowtest/geant/iqmd_Au_Au_1250MeV_1000evts

par1=/cvmfs/hades.gsi.de/install/5.34.34/hydra2-4.9p/defall.sh  # optional par1 : environment script
par2=/lustre/hades/user/sudol/hgeant2-4.9p/hgeant   # optional par2 : executable
par3=""                                                        # optional par3 : geant inidatfile
par4="no"                                                      # optional par4 : 
par5="no"                                                      # optional par5 : number of events
par6="no"                                            # optional par6 "no" -> no vertexfile   "extract_vertex" -> vertex file in outputdir/vertex
par7="no"                                                      # optional par7
resources="--mem=2000 --time=0-8:00:00"                            # runtime < 10h, mem < 2GB

jobarrayFile="gen8_${part}_jobarray.dat"

filelist=${currentDir}/${impact}_new.list  # file list in local dir! not in submissiondir!!!
#filelist=${currentDir}/all_files_${part}_${nFiles}.list  # file list in local dir! not in submissiondir!!!
#filelist=${currentDir}/missing.list  # file list in local dir! not in submissiondir!!!

createList=yes  # (yes/no) use this to create files list with generic names (for simulation, testing)
                 # use "no" if you have a filelist available

#nFiles=$(cat $filelist | wc -l)
######################################################################


#---------------------------------------------------------------------
# create a file list for submission (simulation, testing etc.)
# for real data you will have a filelist with real filenames
if [ "$createList" == "yes" ]
then
   if [ -f $filelist ]
   then
     echo "===> REMOVING EXISTING FILELIST : $filelist"
     rm -f $filelist 
   fi

   echo "===> CREATE FILELIST : $filelist"
   for ((ct=1+$offset;ct<=+$offset+$nFiles;ct++))
   do
      echo ./input/au123au_${part}_${ct}.dat >> $filelist
   done
fi
#---------------------------------------------------------------------

mkdir -p $submissiondir/input

#---------------------------------------------------------------------
# create input inidat for geant
echo $replace -t $template -d ${submissiondir}/input/au123au_${part}_XXX.dat -n $nevents -i ${input1}_XXX.evt -p ${submissiondir}/geom_rich -o ${output}_XXX.root > replaceMask.txt  # pluto lep only
#echo $replace -t $template -d ${submissiondir}/input/au123au_${part}_XXX.dat -n $nevents -i ${input3}_XXX.evt -p ${submissiondir}/geom_rich -o ${output}_XXX_.root > replaceMask.txt  # urqmd only
#echo $replace -t $template -d ${submissiondir}/input/au123au_${part}_XXX.dat -n $nevents -i ${input1}_XXX.evt -u ${input3}_XXX.evt -p ${submissiondir}/geom_rich -o ${output}_XXX_.root > replaceMask.txt  # pluto + urqmd
#echo $replace -t $template -d ${submissiondir}/input/au123au_${part}_XXX.dat -n $nevents -i ${input1}_XXX.evt -s ${input2}_XXX.evt -u ${input3}_XXX.evt -p ${submissiondir}/geom_rich -o ${output}_XXX_.root > replaceMask.txt  # pluto + pluto + urqmd
#---------------------------------------------------------------------

nFiles=$( cat $filelist | wc -l)

#---------------------------------------------------------------------
# create needed dirs
if [ ! -d $submmissionbase ]
then
    echo "===> CREATE SUBMISSIONBASEDIR : $submmissionbase"
    mkdir -p $submmissionbase
else
    echo "===> USE SUBMISSIONBASEDIR : $submmissionbase"
fi

#---------------------------------------------------------------------
# output dirs

if [ ! -d $pathoutputlog ]
then
   echo "===> CREATE LOGDIR : $pathoutputlog"
   mkdir -p $pathoutputlog
else
   echo "===> USE LOGDIR : $pathoutputlog"
fi

if [ ! -d ${outputdir}/vertex ]
then
   echo "===> CREATE VERTEXDIR : ${outputdir}/vertex"
   mkdir -p ${outputdir}/vertex
else
   echo "===> USE VERTEXDIR : ${outputdir}/vertex"
fi

if [ ! -d ${outputdir}/crash ]
then
   echo "===> CREATE CRASHDIR : ${outputdir}/crash"
   mkdir -p ${outputdir}/crash
else
   echo "===> USE CRASHDIR : ${outputdir}/crash"
fi
#---------------------------------------------------------------------




ctF=0          # counter for file number
ctJ=0          # counter for job number
partNumber=0   # counter for part number


if [ -f $jobarrayFile ]
then
  rm -f $jobarrayFile
fi

echo "===> CREATING JOB ARRAY FILE"
#---------------------------------------------------------------------
# read the files list into an job array
declare -a jobarray
ct1=0
for file in $(cat $filelist)
do
   jobarray[$ct1]=$file
   ((ct1+=1))
done
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# loop over the job array and submit parts with
# nFilesPerJob to GE

while ((ctF<$nFiles))
do

     #---------------------------------------------------------------------
     # build comma separated file list
     # per job
     if [ $nFilesPerJob -gt 1 ]
     then
        infileList=${jobarray[${ctF}]}
        ((ctF+=1))
        for (( ctList=1;ctList<$nFilesPerJob; ctList++ ))
        do   	
            if [ $ctF -lt ${nFiles} ]
            then
               infileList="${infileList},${jobarray[${ctF}]}"
               ((ctF+=1))
            fi
        done
     else 
        infileList=${jobarray[${ctF}]}
        ((ctF+=1))
     fi
     #---------------------------------------------------------------------
     
     ((partNumber+=1))

     logfile="${pathoutputlog}/${filename}_${partNumber}.log"

     if [ $nFilesPerJob -eq 1 ]
     then
        file=$(basename ${infileList})
        logfile="${pathoutputlog}/${file}.log"
     fi
     
     if [ -f ${logfile} ]
     then
        rm -f ${logfile}
     fi      
     
     #echo "-----------------------------------------------------------------------------"
     #echo "add part ${partNumber}  last file ${ctF} of $nFiles ====> add new job ${infileList}"

     ######################################################################
     #  SEND NEW JOB (USER SPECIFIC)
     
     par3=${infileList}

           #defall.sh prog filelist outdir  nev   extract_vert
     echo "${par1} ${par2} ${par3} ${par4} ${par5} ${par6} ${par7}" >>  $jobarrayFile
     

     ######################################################################
     
done
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# sync the local modified stuff 
# to the submission dir
echo "===> SYNC CURENTDIR TO SUBMISSIONDIR : rsync  -vHa $currentDir ${submmissionbase}"
rsync  -vHa $currentDir ${submmissionbase}/

syncStat=$?

if [ ! $syncStat -eq 0 ]
then
     echo "===> ERROR : SYNCHRONIZATION ENCOUNTERED PROBLEMS"
else

  echo "-------------------------------------------------"


  nFiles=$( cat $jobarrayFile | wc -l)
  arrayoffset=0;
  ctsend=0
  block=500
  while ((${ctsend} * ${block} < ${nFiles}))
  do
     ((start=${ctsend}*${block}))
     ((rest=${nFiles}-${start}))
     if [ $rest -le $block ]
     then
        ((stop=$rest))
     else
        ((stop=$block))
     fi
     ((arrayoffset=${ctsend} * ${block}))
     command="--array=1-${stop} ${resources} -D ${submissiondir}  --output=${pathoutputlog}/slurm-%A_%a.out ${jobscript} ${submissiondir}/${jobarrayFile} ${pathoutputlog} ${arrayoffset}"
     #echo $command
      sbatch $command

     ((ctsend+=1))
  done

  echo "${nFiles} jobs for part ${part} submitted"
fi
