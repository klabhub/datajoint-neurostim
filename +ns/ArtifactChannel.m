%{
# Trials and time periods that should not be used for analysis, per channel
-> ns.Artifact
channel : int # The channel that has the artifact
---
trial =NULL :      longblob   # A vector of trials that are considered artifacts and should be removed from analysis
start = NULL : longblob  # A vector of neurostim times that indicate the start of an artifact period
stop =NULL  : longblob   # A vector of neurostim times that indicate the end of an artifact period
%}
%
% See also ns.Artifact, which stores the trials and periods that should be
% excluded for all channels. This table only has additional information per
% channel.
classdef ArtifactChannel< dj.Part
     properties (SetAccess = protected)
        master = ns.Artifact
     end 
end