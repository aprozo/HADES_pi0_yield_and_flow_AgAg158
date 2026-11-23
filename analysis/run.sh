# nohup bash run.sh ./parameter_file.conf &

if [ $# -eq 0 ]  # no arguments passed
then
  echo "ERROR: wrong usage of the script. No arguments found"
  echo "       bash run.sh [-r] ./parameter_file.conf "
  echo "       use -r option if you whant to rewrite existing folders"
  exit 15
fi

DO_REMOVE=0
while getopts "r" Option
do
  case $Option in
    r     ) DO_REMOVE=1; echo "    -r option: rewrite existing folders";;
    *     ) echo "    Unknown option.";;   # Default
  esac
done
shift $(($OPTIND - 1))

RUN_URQMD_TO_PT_Y=0

RUN_EXP_DST_TO_PHOTONS=0
RUN_SIM_DST_TO_PHOTONS=0

RUN_PHOTONS_TO_HISTOGRAMS=0
RUN_SIM_FLOW=0
RUN_PHOTONS_TO_HISTOGRAMS_FLOW=0

RUN_HISTOGRAMS_TO_PT_Y=1

list_exp=3200A
#day073
#test
#2500A
#3200A
list_sim=sim
#test
#sim
WHOAMI=`whoami`
if [[ -d /lustre/hades/user ]]
    then
    is_slurm_cluster=1
    else
    is_slurm_cluster=0
fi
#check if all computations finished
wait_slurm()
{
    sleep 20
    if [[ $is_slurm_cluster -eq 1 ]]
        then
        wait_count=0
        wait_step=20
        wait_time=0
        while [[ 1 -lt `squeue -u $WHOAMI | wc -l` && $wait_time -lt 50000 ]]
        do
            echo "Previous task is running for $wait_time seconds. Wait for another 20 sec ($wait_count trial) ..."
            sleep $wait_step
            (( wait_time+=$wait_step ))
            (( wait_count++ ))
        done
        unset wait_count
        unset wait_step
        unset wait_time

        if [[ 1 -ne `squeue -u $WHOAMI | wc -l` ]]
            then
            echo "ABORT: There are unfinished slurm tasks. Cancelling."
            scancel -u $WHOAMI
        fi
    fi
    return 0
}

echo "Start SYSTEMATIC UNCERTAINTY study"
echo ".................................."
echo
parameter_file=$1
if [[ ! -e $parameter_file ]]
    then
    echo "ABORT: file $parameter_file doesn't exist"
    exit 1
fi
parameter_filename=$(basename -- "$parameter_file")
parameter_extension="${parameter_filename##*.}"
if [[ ! "$parameter_extension" = "conf" ]]
    then
    echo "ABORT: wrong extension \"$parameter_extension\" of the parameter file"
    exit 2
fi
#-------------------------------------------------------------------
# Make local output folder
#-------------------------------------------------------------------
parent_dir=`pwd`
output_foldername="${parameter_filename%.*}"
path_local_output_full=$parent_dir/results/$output_foldername
if [[ DO_REMOVE -eq 1 ]]
then
    rm -rf $path_local_output_full
fi
if [ -e "$path_local_output_full" ]
    then
        echo " folder " $path_local_output_full "already exists, rewrite"
    else
        mkdir $path_local_output_full
fi

echo "Output folder: " $path_local_output_full
path_local_sim_to_photons=$path_local_output_full/sim_to_photons
path_local_exp_to_photons=$path_local_output_full/exp_to_photons
#-------------------------------------------------------------------
# Make SLURM directory
#-------------------------------------------------------------------
if [[ $is_slurm_cluster -eq 1 ]]
then
    path_slurm=/lustre/hades/user/${WHOAMI}/mar19/systematic_uncertainty
    echo "Output folder on SLURM cluster: " ${path_slurm}
    if [[ ! -d ${path_slurm} && $is_slurm_cluster -eq 1 ]]
    then
        mkdir ${path_slurm}
        echo "....folder is created"
    fi
    path_slurm_output_full=$path_slurm/$output_foldername
    if [[ RUN_EXP_DST_TO_PHOTONS -ne 0 ]]
    then
        path_slurm_exp_to_photons=$path_slurm_output_full/exp_to_photons
    else
        path_slurm_exp_to_photons=$path_slurm/exp_to_photons
    fi
    if [[ RUN_SIM_DST_TO_PHOTONS -ne 0 ]]
    then
        path_slurm_sim_to_photons=$path_slurm_output_full/sim_to_photons
    else
        path_slurm_sim_to_photons=$path_slurm/sim_to_photons
    fi
    if [[ DO_REMOVE -eq 1 ]]
    then
        rm -rf $path_slurm_output_full
    fi
    if [ -e "$path_slurm_output_full" ]
    then
        echo "folder " $path_slurm_output_full "already exists , rewrite"
    else
        mkdir $path_slurm_output_full
    fi
fi
#-------------------------------------------------------------------
# compile macro parameters_to_csv.C
#-------------------------------------------------------------------
g++ -o parameters_to_csv parameters_to_csv.C `root-config --cflags --libs`
if [[ ! -x parameters_to_csv ]]
    then
    echo "Macro parameters_to_csv.C wasn't properly compiled"
  #  exit 8
fi

# Write all parameters to the summary table
cp $parameter_file $path_local_output_full/parameters.conf
cp GlobalParameters.h $path_local_output_full/

cp $path_local_output_full/parameters.conf $path_slurm_output_full
$parent_dir/parameters_to_csv $parameter_file >> results/summary_table.csv

#===================================================================
# UrQMD to 2D
#===================================================================
path_urqmd_to_pt_y=$path_local_output_full/urqmd_to_pt_y

if [[ $RUN_URQMD_TO_PT_Y -ne 0 ]]
then
    mkdir -p $path_urqmd_to_pt_y
    cp -r $parent_dir/macros/urqmd_to_pt_y/* $path_urqmd_to_pt_y
    cd $path_urqmd_to_pt_y
    ###########
    bash run.sh
    ###########
    if [[ 0 -ne $? ]]
        then
        echo "Making urqmd_pt_y.root failed"
        exit 9
    fi
    cd $parent_dir
else
    echo "Skip urqmd_to_pt_y"
fi


#===================================================================
# DST to photons.root
#===================================================================
#EXPERIMENT
if [[ $RUN_EXP_DST_TO_PHOTONS -ne 0 ]]
then
    while IFS= read -r line
    do
        if [[ $line == "#"* ]]
        then
            continue
        fi
        if [[ $line == *"PATH_TO_DST_EXP"* ]]
        then
        path_to_dst_exp=${line:16}
        fi
    done < $parameter_file
    unset line
    echo "Path to DST files with experimental data: " $path_to_dst_exp
    if [[ ! -d $path_to_dst_exp && $is_slurm_cluster -eq 1 ]]
    then
        echo "ABORT: folder " $path_to_dst_exp " doesn't exist"
        exit 10
    fi
    mkdir -p $path_local_exp_to_photons
    cp -r $parent_dir/macros/exp_sim_to_photons/* $path_local_exp_to_photons
    cd $path_local_exp_to_photons
    if [[ $is_slurm_cluster -eq 1 ]]
    then
        find $path_to_dst_exp/3200A/0{68..88}/root | grep accepted.root > 3200A.list
      #  if [$list_exp=="day073"]
      #  then
      #     find $path_to_dst_exp/3200A/073/root | grep accepted.root > 3200A.list
      #ssh v    fi
        #find $path_to_dst_exp/2500A/???/root | grep accepted.root > 2500A.list
        #find $path_to_dst_exp/200A/???/root | grep accepted.root > 200A.list
        #find $path_to_dst_exp/0A/???/root | grep accepted.root > 0A.list

        cat 3200A.list > full.list
        #cat 200A.list >> full.list
        #cat 0A.list   >> full.list
        head full.list > test.list
    fi
    #########################################################
    # 1 stands for experiment
    # 0 stands for simulation
    bash run.sh 1 $path_slurm_output_full $path_slurm_exp_to_photons $list_exp
    #########################################################
    if [[ 0 -ne $? ]]
    then
        echo "Making photons_exp.root failed"
        exit 11
    fi
    cd $parent_dir
else
    echo "Skip exp_to_photons"
fi


#-------------------------------------------------------------------
#SIMULATION
#-------------------------------------------------------------------

if [[ $RUN_SIM_DST_TO_PHOTONS -ne 0 ]]
then
    wait_slurm
    while IFS= read -r line
    do
        if [[ $line == "#"* ]]
        then
            continue
        fi
        if [[ $line == *"PATH_TO_DST_SIM"* ]]
        then
        path_to_dst_sim=${line:16}
        fi
    done < $parameter_file
    unset line
    echo "Path to DST files with simulation data: " $path_to_dst_sim
    if [[ ! -d $path_to_dst_sim && $is_slurm_cluster -eq 1 ]]
        then
        echo "ABORT: folder " $path_to_dst_sim " doesn't exist"
        exit 12
    fi
    mkdir -p $path_local_sim_to_photons
    cp -r $parent_dir/macros/exp_sim_to_photons/* $path_local_sim_to_photons
    cd $path_local_sim_to_photons
    if [[ $is_slurm_cluster -eq 1 ]]
    then
        find $path_to_dst_sim | grep 19.root > sim.list
        head sim.list > test.list
    fi
    #########################################################
    # 1 stands for experiment
    # 0 stands for simulation
    bash run.sh 0 $path_slurm_output_full $path_slurm_sim_to_photons $list_sim
    #########################################################
    if [[ 0 -ne $? ]]
    then
        echo "Making photons_sim.root failed"
        exit 13
    fi
    cd $parent_dir
else
    echo "Skip sim_to_photons"
fi


#===================================================================
# photons_to_histograms :
#-------------------------------------------------------------------
# Select photon candidates,
# Combine photons to diphotons,
# Fill the histograms
#===================================================================

path_local_photons_to_histograms=${path_local_output_full}/photons_to_histograms
path_slurm_photons_to_histograms=${path_slurm_output_full}/photons_to_histograms

if [[ $RUN_PHOTONS_TO_HISTOGRAMS -ne 0 ]]
then
    wait_slurm
    mkdir -p $path_local_photons_to_histograms
    cp -r $parent_dir/macros/photons_to_histograms/* $path_local_photons_to_histograms
    cd $path_local_photons_to_histograms
    if [[ $is_slurm_cluster -eq 1 ]]
    then
        if [[ $RUN_EXP_DST_TO_PHOTONS -ne 0 ]]
        then
            find $path_slurm_exp_to_photons/$list_exp | sort -n | grep root > ${list_exp}_exp.list
        else
            find /lustre/hades/user/${WHOAMI}/mar19/systematic_uncertainty/exp_to_photons/$list_exp | sort -n | grep root > ${list_exp}_exp.list
        fi
        head ${list_exp}_exp.list > test_exp.list
        if [[ $RUN_SIM_DST_TO_PHOTONS -ne 0 ]]
        then
            find $path_slurm_sim_to_photons/$list_sim | sort -n |grep root > ${list_sim}_sim.list
        else
            find /lustre/hades/user/${WHOAMI}/mar19/systematic_uncertainty/sim_to_photons/$list_sim | sort -n | grep root > ${list_sim}_sim.list
        fi
        head ${list_sim}_sim.list > test_sim.list
    fi

    # run macros
    ###################################################################################
    # experiment
    name_photons_to_histograms_exp=${list_exp}_exp
    echo "bash run.sh 1 $path_slurm_output_full $path_slurm_photons_to_histograms $name_photons_to_histograms_exp" > command.txt
    bash run.sh 1 $path_slurm_output_full $path_slurm_photons_to_histograms $name_photons_to_histograms_exp
    status1=$?
    ###################################################################################
    # simulation
    wait_slurm
    name_photons_to_histograms_sim=${list_sim}_sim
    echo "bash run.sh 0 $path_slurm_output_full $path_slurm_photons_to_histograms $name_photons_to_histograms_sim" >> command.txt
    bash run.sh 0 $path_slurm_output_full $path_slurm_photons_to_histograms $name_photons_to_histograms_sim
    status2=$?
    ###################################################################################
    wait_slurm
    # check if the scripts are successfully run
    if [[ 0 -ne status1 || 0 -ne status2 ]]
    then
        echo "Combination of diphotons failed"
  #      exit 14
    fi

    if [[ $is_slurm_cluster -eq 1 ]]
    then
        cd $path_local_photons_to_histograms
        hadd -j -f -k ${list_exp}_exp.root `find $path_slurm_photons_to_histograms/${list_exp}_exp | grep root`
        hadd -j -f -k ${list_sim}_sim.root `find $path_slurm_photons_to_histograms/${list_sim}_sim | grep root`
        cd $parent_dir
    fi
else
    echo "Skip photons_to_diphotons"
fi






#===================================================================
# photons_to_histograms :
#-------------------------------------------------------------------
# Select photon candidates,
# Combine photons to diphotons,
# Fill the histograms
#===================================================================

path_local_photons_to_histograms_flow=${path_local_output_full}/photons_to_histograms_flow
path_slurm_photons_to_histograms_flow=${path_slurm_output_full}/photons_to_histograms_flow

if [[ $RUN_PHOTONS_TO_HISTOGRAMS_FLOW -ne 0 ]]
then
    wait_slurm
    mkdir -p $path_local_photons_to_histograms_flow
    cp -r $parent_dir/macros/photons_to_histograms_flow/* $path_local_photons_to_histograms_flow
    cd $path_local_photons_to_histograms_flow
    if [[ $is_slurm_cluster -eq 1 ]]
    then
        if [[ $RUN_EXP_DST_TO_PHOTONS -ne 0 ]]
        then
            find $path_slurm_exp_to_photons/$list_exp | sort -n | grep root > ${list_exp}_exp.list
        else
            find /lustre/hades/user/${WHOAMI}/mar19/systematic_uncertainty/exp_to_photons/$list_exp | sort -n | grep root > ${list_exp}_exp.list
        fi
        head ${list_exp}_exp.list > test_exp.list
        if [[ $RUN_SIM_DST_TO_PHOTONS -ne 0 ]]
        then
            find $path_slurm_sim_to_photons/$list_sim | sort -n |grep root > ${list_sim}_sim.list
        else
            find /lustre/hades/user/${WHOAMI}/mar19/systematic_uncertainty/sim_to_photons/$list_sim | sort -n | grep root > ${list_sim}_sim.list
        fi
        head ${list_sim}_sim.list > test_sim.list
    fi

    # run macros
    ###################################################################################
    # experiment
    name_photons_to_histograms_flow_exp=${list_exp}_exp
    echo "bash run.sh 1 $path_slurm_output_full $path_slurm_photons_to_histograms_flow $name_photons_to_histograms_flow_exp" > command.txt
    bash run.sh 1 $path_slurm_output_full $path_slurm_photons_to_histograms_flow $name_photons_to_histograms_flow_exp
    status1=$?

    # check if the scripts are successfully run
    if [[ 0 -ne status1 ]]
    then
        echo "Combination of diphotons failed"
    fi

    if [[ $is_slurm_cluster -eq 1 ]]
    then
        wait_slurm
        cd $path_local_photons_to_histograms_flow
        hadd -j -f -k ${list_exp}_exp.root `find $path_slurm_photons_to_histograms_flow/${list_exp}_exp | grep root`
        cd $parent_dir
    fi

    ###################################################################################

    # simulation
    if [[ $RUN_SIM_FLOW -ne 0 ]]
    then
    name_photons_to_histograms_flow_sim=${list_sim}_sim
    echo "bash run.sh 0 $path_slurm_output_full $path_slurm_photons_to_histograms_flow $name_photons_to_histograms_flow_sim" >> command.txt
    bash run.sh 0 $path_slurm_output_full $path_slurm_photons_to_histograms_flow $name_photons_to_histograms_flow_sim
    ###################################################################################
    wait_slurm

     if [[ $is_slurm_cluster -eq 1 ]]
     then
        cd $path_local_photons_to_histograms_flow
        hadd -j -f -k ${list_sim}_sim.root `find $path_slurm_photons_to_histograms_flow/${list_sim}_sim | grep root`
        cd $parent_dir
     fi

    fi


else
    echo "Skip photons_to_diphotons"
fi


#===================================================================
# histograms_to_pt_y :
#-------------------------------------------------------------------
# Extract pi0 number from diphoton histograms
#===================================================================

path_local_histograms_to_pt_y=${path_local_output_full}/histograms_to_pt_y
path_slurm_histograms_to_pt_y=${path_slurm_output_full}/histograms_to_pt_y

if [[ $RUN_HISTOGRAMS_TO_PT_Y -ne 0 ]]
then
    mkdir -p $path_local_histograms_to_pt_y
    cp -r $parent_dir/macros/histograms_to_pt_y/* $path_local_histograms_to_pt_y
    cd $path_local_histograms_to_pt_y

    # run macros
    ###################################################################################
    input_hists_sim=${path_local_photons_to_histograms}/${list_sim}_sim.root
    input_hists_exp=${path_local_photons_to_histograms}/${list_exp}_exp.root



    if [ ! -f "$input_hists_exp" ]; then
        echo "exp hists do not exist"
        exit 15
    fi

    if [ ! -f "$input_hists_sim" ]; then
        echo "sim hists do not exist"
        exit 16
    fi
    if [ ! -f "${path_urqmd_to_pt_y}/urqmd_pt_y.root" ]; then
        echo "urqmd hists do not exist"
        exit 17
    fi

    root -l -n -b -q 'histograms_to_pt_y.cpp++("'$input_hists_exp'","'$input_hists_sim'","'${path_urqmd_to_pt_y}'/urqmd_pt_y.root")'

    #if [ $RUN_PHOTONS_TO_HISTOGRAMS_FLOW -ne 0 ];
    #then
      input_hists_sim_flow=${path_local_photons_to_histograms_flow}/${list_sim}_sim.root
      input_hists_exp_flow=${path_local_photons_to_histograms_flow}/${list_exp}_exp.root
      root -l -n -b -q 'histograms_to_pt_y_flow.cpp++("../photons_to_histograms_flow/3200A_exp.root")'
    #fi
    echo "end of run for paramters set - $parameter_file"
###################################################################################
else
    echo "Skip histograms_to_pt_y"
fi

# Final results
echo
echo ".................................."
echo "Stop  SYSTEMATIC UNCERTAINTY study"
exit 0
