%{
# A mapping that determines which experiments belong to a Preprocessed set.
-> sbx.Preprocessed     # Preprocessed Set
-> ns.Experiment
%}

classdef PreprocessedExperiment< dj.Part
    properties (SetAccess = protected)
        master = sbx.Preprocessed
    end

end
