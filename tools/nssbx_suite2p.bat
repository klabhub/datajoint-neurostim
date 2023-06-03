@echo off
REM  Windows Batch File Wrapper to run suite2p 
REM This is called from sbx.Preprocessed.makeTuples with 5 input arguments:
REM 1 - The miniconda activation batch file
REM 2 - The miniconda install folder
REM 3 - The python wrapper module (nssbx_suite2p.py) 
REM 4 - The temporary file with the ops input for suite2p
REM 5 - The temporary file with the db input for suite2p

REM Activate Conda Suite2p, then call python wrapper.
CALL %1 %2  &  conda activate suite2p & python %3 %4 %5