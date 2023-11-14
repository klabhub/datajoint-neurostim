#!/bin/bash 
# Bash script Wrapper to run suite2p 
# This is called from sbx.Preprocessed.makeTuples with 3 input arguments:
# 1 - The installation location of conda (e.g. ~/miniconda3)
# 2 - The python wrapper module (nssbx_suite2p.py) 
# 3 - The temporary file with the ops input for suite2p
# 4 - The temporary file with the db input for suite2p


# Because this runs in a fresh bash shell, we need to init the conda hook
eval "$(conda shell.bash hook)"
source "$1/etc/profile.d/conda.sh"

# Now we can activate the suite2p environment
conda activate suite2p
# Run the python script
python $2 $3 $4