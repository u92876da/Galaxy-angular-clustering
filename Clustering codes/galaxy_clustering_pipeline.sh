#!/bin/bash

# Galaxy clustering pipeline .sh
# calculate things if not already done so in this script
# this script is for power law fits only

# ENVIRONMENT variables ----------------------------------
CAT_FILENAME="catalogue"
# 
# what to plot
PLOT_W_THETA=true
PLOT_BIAS=true
# etc to plot all sorts of things


export CAT_FILENAME
export PLOT_W_THETA
export PLOT_BIAS

# a big if statement in here to plot the relevant things
python galaxy_field.py


# 1) calculate w_theta ACFs or CCFs for relevant galaxy populations and save the data
# 2) fit to a power law and calculate biases and save the parameters


# store data for each individually as well as combined into presentable graphs/tables etc
