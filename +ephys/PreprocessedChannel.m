%{
# Preprocessed signals per channel 
-> ephys.Preprocessed
channel : int # The channel that recorded the signal
---
signal : longblob # The preprocessed signal for each sample in the experiment
info = NULL : longblob # Any additional information on this preprocessed channel
%}
%
% See also ephys.Preprocessed
classdef PreprocessedChannel < dj.Part
     properties (SetAccess = protected)
        master = ephys.Preprocessed;  
    end

    methods  (Access = protected)
        function makeTuples(~,~)
             %Handled by the parent class ephys.Preprocessed
        end
    end
end