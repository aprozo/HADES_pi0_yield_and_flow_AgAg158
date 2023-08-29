#------------------------
# runs macro which draws
# 2D pt-y histogram from
# UrQMD data
#------------------------
echo
while IFS= read -r line
do
    if [[ $line == "#"* ]]
    then
        continue
    fi
    if [[ $line == *"PATH_TO_CORRECTION_DATA"* ]]
    then
        path_to_correction_data=${line:24}
    fi
    if [[ $line == *"PATH_TO_URQMD_DATA"* ]]
    then
        path_to_urqmd_data=${line:19}
    fi
done < ../parameters.conf
unset line

echo "Start convert UrQMD data to 2D Pt - Y plot"
echo ".........................................."
if [ -e "urqmd_to_pt_y" ]
then
    rm urqmd_to_pt_y
fi
if [ ! -e "$path_to_urqmd_data" ]
then
    echo "ERROR: no file $path_to_urqmd_data"
    exit 2
else
    echo "UrQMD data: $path_to_urqmd_data"
fi
if [ ! -e "$path_to_correction_data/acceptancePhiThetaMar19.root" ]
then
    echo "ERROR: no file $path_to_correction_data/acceptancePhiThetaMar19.root"
    exit 2
else
    echo "Acceptance of Ecal: $path_to_correction_data/acceptancePhiThetaMar19.root"
fi
g++ -o urqmd_to_pt_y urqmd_to_pt_y.cc `root-config --cflags --libs`
if [ ! -e "urqmd_to_pt_y" ]
    then
    echo "ERROR: urqmd_to_pt_y.cc wasn't properly compiled"
    #exit 3
fi

./urqmd_to_pt_y 100000000
exit_code=$?
if [[ $exit_code -ne 0 ]]
then
    echo "ERROR: UrQMD wasnt properly converted to 2D histograms"
    echo "Exit code: $exit_code"
fi
echo ".........................................."
echo "Stop  convert UrQMD data to 2D Pt - Y plot"
echo
#exit $exit_code
