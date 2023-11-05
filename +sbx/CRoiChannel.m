%{
# Preprocessed continuous signals per channel 
-> sbx.CRoi
channel : int # The channel that recorded the signal
---
signal : longblob # The preprocessed signal for each sample in the experiment
name = NULL : varchar(128) # A name for this channel 
info = NULL : longblob # Information on this preprocessed channel
%}
%
% spikes, fluorescence  neuropil
% See also sbx.CRoi
classdef CRoiChannel < dj.Part
     properties (SetAccess = protected)
        master = sbx.CRoi
     end 
end