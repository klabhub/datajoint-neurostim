#!/bin/bash 
# Bash script Wrapper to run NWB/Dandi 
# This is called from ns/Experiment/nwbExport with input arguments:
# 1 - The installation location of conda (e.g. ~/miniconda3)
# 2 - The conda environment for NWB/Dandi 
# 3 - The command to run 


__conda_setup="$('$1/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "$1/etc/profile.d/conda.sh" ]; then
        source "$1/etc/profile.d/conda.sh"
    else
       export PATH="$1/bin:$PATH"
    fi
fi
unset __conda_setup


# Now we can activate the environment
conda activate $2
# Run the python script
eval "$3"
