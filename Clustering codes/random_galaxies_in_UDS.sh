#!/bin/bash


# random galaxies in UDS .sh


# in .py file; ngals = os.environ["NGALS"]

start=$SECOND

# while getopts i:o: flag
# do
#      case "${flag}" in
#          i) MASK_DIR=${OPTARG};
#          o) SAVE_DIR=${OPTARG};
#      esac
#  done

NGALS=1000
RANDOM_GALS_FILE="UDS_random_gals"
MASK_DIR=$1

export NGALS
export RANDOM_GALS_FILE
export MASK_DIR

python two_pt_angular_correlation.py

# e.g. pip install astropy 

duration=$((SECONDS - start))
echo duration

# for loop requires starting and closing with: do; done