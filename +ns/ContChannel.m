%{
# Preprocessed continuous signals per channel 
-> ns.Cont
channel : int # The channel that recorded the signal
---
signal : longblob # The preprocessed signal for each sample in the experiment
name = NULL : varchar(128) # A name for this channel 
info = NULL : longblob # Information on this preprocessed channel
%}
%
% See also ns.Cont
classdef ContChannel < dj.Part
     properties (SetAccess = protected)
        master = ns.Cont  
     end 
end