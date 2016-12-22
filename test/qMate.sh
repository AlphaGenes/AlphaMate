#!/bin/bash

########################################
# GE job script for ECDF Cluster       #
# ecdf-systems-team@lists.ed.ac.uk     #
########################################

# Grid Engine options
#$ -cwd

# Initialise the module environment
source /etc/profile.d/modules.sh

# Load modules
module add R # for cpumemlogplot.R
module add intel/2016

# Standard report
echo
pwd
date

echo
echo "System PATH (default):"
echo $PATH
echo "System PATH modified:"
export PATH=.:~/bin:$PATH
echo $PATH
echo "Ulimit:"
ulimit -a
echo

echo "Starting job:"
time ../bin/AlphaMate &

# CPU and RAM tracking                                                                 
cpumemlog.sh $!

echo
pwd
date

