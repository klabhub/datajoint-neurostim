%{
# Preprocessed and epoched data
-> ns.Epoch
channel : int
---
signal : longblob
%}

classdef EpochChannel < dj.Part

    properties

        master = ns.Epoch

    end
    
end