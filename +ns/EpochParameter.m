%{

# Epoching parameters to segment trial data. All times are relative to the plugin time

etag : varchar (32) # unique tag
paradigm : varchar(32)
dimension : varchar(32) # Condition from the dimension table
---
plugin_name : varchar(32) # The plugin to be timelocked to from DimensionCondition table
pv : BLOB # structure array containing the custom parameters

%}
%
% parms
%   parms.epoch_win double (1,2) # The relative start and end time of the epochs with respect to the plugin time
%   parms.baseline bool (1,1) # whether to baseline or not
%   parms.baseline_win double (1,2) # The baseline window (1, 1)
%
% MOz Feb, 2025
classdef EpochParameter < dj.Lookup

end