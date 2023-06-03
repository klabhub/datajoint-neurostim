#!/bin/bash 
# Bash script Wrapper to run suite2p 
# This is called from sbx.Preprocessed.makeTuples with 3 input arguments:
# 1 - The python wrapper module (nssbx_suite2p.py) 
# 2 - The temporary file with the ops input for suite2p
# 3 - The temporary file with the db input for suite2p

# activate Conda Suite2p
conda activate suite2p 

python $1 $2 $3