%{
#  Preprocessing instructions used for continuous data
tag         :  varchar(128)     # A  unique name for these preprocessing instructions
---
fun         : varchar(255)      # The function that reads the data and does the preprocessing. 
extension   : varchar(10)       # Extension of the data files to which this applies. 
description : varchar(1024)     # Short description
parms       : longblob          # struct containing all parameters that the fun needs to read the file and preprocess the data.
include=NULL : varchar(128)      # wilcard for filenames to include
exclude=NULL : varchar(128)      # wildcard for filenames to exclude
%}
%
% EXAMPLE
% struct('tag','eeg','fun','ephys.intan.read',
%           'description','EEG  preprocessing for Intan RHS files',
%           'extension','.rhs','parms',struct('downsample',1000,'channel',1:32))
% See ephys.intan.read to understand what 'parms' can contain.
%
% BK - June 2023

classdef ContParm < dj.Lookup

end