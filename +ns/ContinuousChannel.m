%{
# Preprocessed continuous signals per channel 
-> ns.Continuous
channel : int # The channel that recorded the signal
---
signal : longblob # The preprocessed signal for each sample in the experiment
info = NULL : longblob # Any additional information on this preprocessed channel
%}
%
% See also ns.Continuous
classdef ContinuousChannel < dj.Part
     properties (SetAccess = protected)
        master = ns.Continuous;  
     end 
end