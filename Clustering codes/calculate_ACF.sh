#!/bin/bash
export FROM_SHELL=true

# calculate w_theta ACF for specific galaxy populations

export FIELDNAME="Rachana UDS DR11"
Z_MID_ARR=(1 2 3 4)
export Z_WIDTH=0.6
STELLAR_MASS_MIN_ARR=("10" "10.5" "") # in units of log_10
STELLAR_MASS_MAX_ARR=("" "" "") # in units of log_10
#M_UV_MAX=("") # i.e. dimmest things
export CALCULATE_ACF=true
echo $CALCULATE_ACF

cd /Users/user/Documents/PGR
mkdir -p "$FIELDNAME/w_theta ACF/Raw data"

z_length=${#Z_MID_ARR[@]}
stellar_mass_length=${#STELLAR_MASS_MIN_ARR[@]}
for ((i=0; i<${z_length}; i++));
do
    cd "/Users/user/Documents/PGR/$FIELDNAME/w_theta ACF/Raw data"
    export Z_MID=${Z_MID_ARR[$i]}
    mkdir -p "z=$Z_MID, dz=$Z_WIDTH"

    for ((j=0; j<${stellar_mass_length}; j++));
    do
        cd "/Users/user/Documents/PGR/$FIELDNAME/w_theta ACF/Raw data/z=$Z_MID, dz=$Z_WIDTH"
        export STELLAR_MASS_MIN=${STELLAR_MASS_MIN_ARR[$j]}
        export STELLAR_MASS_MAX=${STELLAR_MASS_MIN_ARR[$j]}
        #echo $STELLAR_MASS_MIN
        
        python "/Users/user/Documents/PGR/Clustering codes/two_pt_angular_correlation.py" # calculate ACF
        
    done
done
    



