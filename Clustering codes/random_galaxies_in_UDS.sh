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
RANDOM_GALS_FILE="UDS random gals ra dec"
RANDOM_GALS_DIR="/Users/user/Documents/PGR/UDS field"
CODE_DIR="/Users/user/Documents/PGR/Galaxy clustering/Clustering codes"
MASK_DIR="/Users/user/Documents/PGR/UDS field"

export NGALS
export RANDOM_GALS_FILE
export MASK_DIR

cd "$RANDOM_GALS_DIR"
python "$CODE_DIR/two_pt_angular_correlation.py"

# e.g. pip install astropy 

duration=$((SECONDS - start))
echo duration

# for loop requires starting and closing with: do; done