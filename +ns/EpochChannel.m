%{
# Preprocessed and epoched data
-> ns.Epoch
channel : int
---
signal : longblob
%}

classdef EpochChannel < dj.Part

    properties (SetAccess = protected)

        master = ns.Epoch

    end
    
end