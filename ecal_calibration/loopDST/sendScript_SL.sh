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
#                                                      
# job log files will be named like inputfiles. If nFilesPerJob > 1
# the log files will contain the partnumber.
#
######################################################################
#   CONFIGURATION

user=$(whoami)
currentDir=$(pwd)
day=${1}
acc_or_rej=${2}
sys=${3}

submmissionbase=/lustre/hades/user/${user}/sub/mar19/gen5/calibration
submissiondir=${submmissionbase}/loopDST
nFilesPerJob=50 # for SIM or EXP extract
jobscript=${submissiondir}/jobScript_SL.sh                        # exec script (full path, call without dot, set it executable!)
outputdir=/lustre/hades/user/${user}/mar19/leptonttree/gen5/${sys}_0${day} #${day}   # outputdir for files AND logFiles
pathoutputlog=${outputdir}/out                                  # protocol from batch farm for each file
filename=${acc_or_rej}
jobarrayFile="jobarray_${sys}_0${day}_${acc_or_rej}.dat"
filelist=lists/day_0${day}_${sys}_${acc_or_rej}.list

#
# nFilesPerJob=5 # for SIM or EXP extract
# jobscript=${submissiondir}/jobScript_SL.sh                        # exec script (full path, call without dot, set it executable!)
# outputdir=/lustre/hades/user/${user}/mar19/ChargedParticlesTime/   # outputdir for files AND logFiles
# pathoutputlog=${outputdir}/out                                  # protocol from batch farm for each file
# filename=chargedParticles
# jobarrayFile="jobarray_chargedParticles.dat"
# filelist=lists/day_080_ag158ag_3200A_accepted.list








par1=/cvmfs/hadessoft.gsi.de/install/debian12/6.36.02/hydra2-7.3/defall.sh  # optional par1 : environment script
export SLURM_SINGULARITY_CONTAINER=/cvmfs/vae.gsi.de/vae25/containers/user_container-production.sif  # Debian 12 container for debian12 hydra
par2=${submissiondir}/loopDST                          # optional par2 : executable
par3=""                                                         # optional par3 : input file list
par4=${outputdir}                                               # optional par4 : outputfile (part number will be appended (_num.root))
par5=100000000                                                  # optional par5 : number of events
par6=                                                           # optional par6
par7=                                                           # optional par7
resources="--mem=2500 --time=0-4:00:00"                         # runtime < 10h, mem < 2GB

pathoutputlog=${outputdir}/out

######################################################################



createJobarray=yes # create jobarry if yes
#---------------------------------------------------------------------
# create a file list for submission (simulation, testing etc.)
# for real data you will have a filelist with real filenames
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

if [ ! -d $outputdir ]
then
   echo "===> CREATE OUTPUTDIR : $outputdir"
   mkdir -p $outputdir
else
   echo "===> USE OUTPUTDIR : $outputdir"
fi

#if [ ! -d $outputdir/crash ]
#then
#   echo "===> CREATE CRASHDIR : $outputdir/crash"
#   mkdir -p $outputdir/crash
#fi
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
        par4=${outputdir}/${filename}_${file}.root
     else
        par4=${outputdir}/${filename}_${partNumber}.root
     fi

     if [ -f ${logfile} ]
     then
        rm -f ${logfile}
     fi

     ######################################################################
     #  SEND NEW JOB (USER SPECIFIC)

     par3=${infileList}

     #     defall.sh prog  filelist outfile  nev
     echo "${par1} ${par2} ${par3} ${par4} ${par5} ${par6} ${par7}" >>  $jobarrayFile


     ######################################################################

done
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# sync the local modified stuff
# to the submission dir
echo "===> SYNC CURENTDIR TO SUBMISSIONDIR : rsync  -vHaz $currentDir ${submmissionbase}"
rsync  -vHaz $currentDir ${submmissionbase}/

syncStat=$?

if [ ! $syncStat -eq 0 ]
then
     echo "===> ERROR : SYNCHRONIZATION ENCOUNTERED PROBLEMS"
else

  echo "-------------------------------------------------"

  nFiles=$( cat $jobarrayFile | wc -l)
  ctsend=0
  block=1000
  while ((${ctsend} * ${block} < ${nFiles}))
  do
     ((start=${ctsend}*${block}+1))
     ((stop= ${start}+${block}-1))
     ((rest=${nFiles}-${start}))
     if [ $rest -le $block ]
     then
        ((stop=$start+$rest))
     fi
     command="--array=1-${stop} ${resources} -D ${submissiondir}  --output=${pathoutputlog}/slurm-%A_%a.out -- ${jobscript} ${submissiondir}/${jobarrayFile} ${pathoutputlog} ${arrayoffset}"
     # new
     #command="--array=1-${stop} ${resources} -D ${submissiondir}  --output=${pathoutputlog}/slurm-%A_%a.out -- ${wrap} ${jobscript} ${submissiondir}/${jobarrayFile} ${pathoutputlog} ${arrayoffset}"
     sbatch $command

     ((ctsend+=1))
  done

  echo "${nFiles} jobs for day ${day} submitted"
fi
