#!/bin/bash

export FROM_SHELL=true
export ORIG_DIR="/Users/user/Documents/PGR/Galaxy clustering"
export CALCULATE_ACF=true
export PLOT_ACF=true
export FIT_POWER_LAW=true
export CALCULATE_BIAS=true # bias, and r_0, σ_8_gal intermediate steps
export PLOT_Z_DISTRIBUTION=true

export FIELDNAME="Rachana UDS DR11"
export FILENAME="/Users/user/Documents/PGR/UDS field/DR11-2arcsec-Jun-30-2019.fits"
export RANDOM_GALS_FILE="UDS random gals ra dec.csv"
Z_MID_ARR=(4)
export Z_WIDTH=0.6
STELLAR_MASS_MIN_ARR=("10") # in units of log_10
STELLAR_MASS_MAX_ARR=("") # in units of log_10
#M_UV_MAX=("") # i.e. dimmest things
#COLOUR=("") # red or blue via UVJ selection
export DELTA=0.8 # set as "" if left as free parameter

cd "$ORIG_DIR"
mkdir -p "$FIELDNAME/w_theta ACF/Raw data"
# mkdir -p "$FIELDNAME/w_theta ACF/Individual plots"
# mkdir -p "$FIELDNAME/w_theta ACF/Power law parameters"
# mkdir -p "$FIELDNAME/w_theta ACF/Redshift distribution plots"


z_length=${#Z_MID_ARR[@]}
stellar_mass_length=${#STELLAR_MASS_MIN_ARR[@]}
for ((i=0; i<${z_length}; i++));
do
    export Z_MID=${Z_MID_ARR[$i]}
    echo "z = $Z_MID"

    cd "$ORIG_DIR/$FIELDNAME/w_theta ACF/Raw data"
    mkdir -p "z=$Z_MID, dz=$Z_WIDTH"
    #cd "/Users/user/Documents/PGR/Galaxy clustering/$FIELDNAME/w_theta ACF/Individual plots"
    #mkdir -p "z=$Z_MID, dz=$Z_WIDTH"
    #cd "/Users/user/Documents/PGR/Galaxy clustering/$FIELDNAME/w_theta ACF/Power law parameters"
    #mkdir -p "z=$Z_MID, dz=$Z_WIDTH"

    for ((j=0; j<${stellar_mass_length}; j++));
    do
        export STELLAR_MASS_MIN=${STELLAR_MASS_MIN_ARR[$j]}
        export STELLAR_MASS_MAX=${STELLAR_MASS_MAX_ARR[$j]}
    
        # calculate ACF
        cd "$ORIG_DIR/$FIELDNAME/w_theta ACF/Raw data/z=$Z_MID, dz=$Z_WIDTH"
        python "$ORIG_DIR/Clustering codes/galaxy_field.py"
        
        # calculate power law parameters
        # TO DO: make clear warning sign to calculate the ACF if it doesn't already exist
        #cd "/Users/user/Documents/PGR/Galaxy clustering/$FIELDNAME/w_theta ACF/Power law parameters/z=$Z_MID, dz=$Z_WIDTH"
        
        
    done
done
    



