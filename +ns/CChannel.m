%{
# Preprocessed continuous signals per channel 
-> ns.C
channel : int # The channel that recorded the signal
---
signal : longblob # The preprocessed signal for each sample in the experiment
name = NULL : varchar(128) # A name for this channel 
channelinfo = NULL : longblob # Information on this preprocessed channel
%}
%
% See also ns.C
classdef CChannel < dj.Part
     properties (SetAccess = protected)
        master = ns.C  
     end 
end