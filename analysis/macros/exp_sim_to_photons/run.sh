#par1: isExperiment (1 or 0)
#par2: slurm/../systematic uncertainty
#par3: exp_to_photons or sim_to_photons
#par4: name of list of input files
#-------------------------------------------------
# Macros which creates photons tree from dst files
#-------------------------------------------------
WHOAMI=`whoami`
submmissionbase=${2}
submissiondir=${3}
list_name=${4}
if [[ $1 -eq 0 ]]
then
    echo "Start convert simulation   DST files to photons.root tree"
    echo "........................................................."
    #sim
    #test
    is_experiment=0
else
    echo "Start convert experimental DST files to photons.root tree"
    echo "........................................................."
    #2500A
    #200A
    #0A
    #full
    #test
    is_experiment=1
fi

currentDir=$(pwd | xargs -i basename {})
jobscript=${submissiondir}/jobScript_SL.sh
outputdir=${submissiondir}/${list_name}
path_log=${outputdir}/log
path_crash=${outputdir}/crash

#par1=/cvmfs/hadessoft.gsi.de/install/debian8/install/5.34.34/hydra2-5.6/defall.sh
#par1=/cvmfs/hades.gsi.de/install/5.34.34/hydra2-5.6/defall.sh
#par1=/cvmfs/hades.gsi.de/install/6.12.06/hydra2-5.6/defall.sh
par1=/cvmfs/hadessoft.gsi.de/install/debian12/6.36.02/hydra2-7.3/defall.sh
export SLURM_SINGULARITY_CONTAINER=/cvmfs/vae.gsi.de/vae25/containers/user_container-production.sif  # Debian 12 container for debian12 hydra
#/cvmfs/hades.gsi.de/install/5.34.34/hydra2-5.6/defall.sh # hydra release
par2=${submissiondir}/exp_sim_to_photons                      # executable
par3=""                                                       # input file (to be done later)
par4=""                                                       # output filename (to be done later)
par5=${is_experiment}                                         # bool parameter
par6=-1                                                       # number of events (-1 to analize all)
par7="no"                                                     #
resources="--mem=8000 --time=0-6:00:00"                       #
jobarrayFile="jobarray.dat"                                   #
filelist=./${list_name}.list                                  #

nFiles=$( cat $filelist | wc -l)
nFilesPerJob=20

#---------------------------------------------------------------------
# create needed dirs
if [ ! -d $submissiondir ]
then
    echo "===> CREATE SUBMISSIONDIR : $submissiondir"
    mkdir -p $submissiondir
else
    echo "===> USE SUBMISSIONDIR : $submissiondir"
fi

if [ ! -d $outputdir ]
then
   echo "===> CREATE OUTPUTDIR : $outputdir"
   mkdir -p $outputdir
else
    echo "ABORT: " $outputdir " already exists"
    exit 1
fi

if [ ! -d $path_log ]
then
   echo "===> CREATE LOGDIR : $path_log"
   mkdir -p $path_log
else
   echo "===> USE LOGDIR : $path_log"
fi

if [ ! -d $path_crash ]
then
   echo "===> CREATE CRASHDIR : $path_crash"
   mkdir -p $path_crash
else
   echo "===> USE CRASHDIR : $path_crash"
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
    #---------------------------------------------------------------------

    ((partNumber+=1))

    logfile="${path_log}/${partNumber}.log"

    par4=${outputdir}/${partNumber}.root

    if [ -f ${logfile} ]
    then
       rm -f ${logfile}
    fi

    ######################################################################
    #  SEND NEW JOB (USER SPECIFIC)
    par3=${infileList}
    #     defall.sh prog infilelist outfile isexperiment nev
    echo "${par1} ${par2} ${par3} ${par4} ${par5} ${par6} ${par7}" >>  $jobarrayFile
    ######################################################################
done
#---------------------------------------------------------------------
make
rm -rf build/
echo "===> SYNC CURENTDIR TO SUBMISSIONDIR : rsync  -vHaz $currentDir ${submmissionbase}"
rsync  -vHaz ../$currentDir ${submmissionbase}/
syncStat=$?
if [ ! $syncStat -eq 0 ]
then
    echo "===> ERROR : SYNCHRONIZATION ENCOUNTERED PROBLEMS"
else
    echo "-------------------------------------------------"
    nFiles=$( cat $jobarrayFile | wc -l)
    arrayoffset=0;
    ctsend=0
    block=1000
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
        command="--array=1-${stop} ${resources} -D ${submissiondir}  --output=${path_log}/slurm-%A_%a.out -- ${jobscript} ${submissiondir}/${jobarrayFile} ${path_log} ${arrayoffset}"
        echo $command
        sbatch $command
        ((ctsend+=1))
    done
    echo "${nFiles} jobs submitted"
fi

sleep 10
squeue -u $WHOAMI

echo "........................................................."
if [[ $1 -eq 0 ]]
then
    echo "Stop  convert simulation   DST files to photons.root tree"
else
    echo "Stop  convert experimental DST files to photons.root tree"
fi
echo

exit 0
