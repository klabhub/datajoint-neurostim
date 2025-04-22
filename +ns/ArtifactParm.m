%{
#  Preprocessing instructions used to detect artifact trials and time periods.
atag         :  varchar(32)     # A  unique name for these instructions
ctag         : varchar(32)    # The ctag of the data to which this applies. 
---
fun         : varchar(255)      # The user defined function that does the work. 
description : varchar(1024)     # Short description
parms       : longblob          # struct containing all parameters that the fun needs to do its job.
%}
%
% EXAMPLE
% struct('atag','outliers','fun','artifacts',
%           'description','Finding outliers in EEG data',
%           'ctag','eeg','parms',struct('z',5,'maxVoltage',50e-6))
% This will call the function artifacts.m for each row in C that has the
% ctag 'eeg'. The artifacts.m function has the following prototype:
% 
% [perExpt,perChannel] = fun(C,parms)
% 
%  When populate(ns.Artifact) is called, the fun will be called with C a
%  relevant row from the ns.C table, and parms the parms defined in this
%  table.
% The user defined function should return two structs (which can be empty)
% perExpt defines trials (.trial) and time periods (.start, .stop) that
% identify trials that shoudl be excluded from analysis (in ns.C/align)
% irrespective of the channel.
% perChannel is a struct array that defines the trials (.trial) and the
% tiem periods (.start , .stop) in which a specific channel (.channel)
% should be excluded from analysis.
%
% For an example, see ephys.artifactDetection.m
% See Also ns.Artifact ns.ArtifactChannel ephys.artifactDetection
% BK - December 2023
classdef ArtifactParm < dj.Lookup

end