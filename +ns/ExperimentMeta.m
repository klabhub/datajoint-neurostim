%{
# Meta data for a Neurostim experiment. 
-> ns.Experiment
meta_name : varchar(255)     # The name of the meta data
---
meta_value : varchar(2048)       # The value of the meta data
%}
% nsScan reads from JSON and nsAddToDataJoint inserts this information into
% this table. Use the nsMeta app to update meta data.
%
% BK - Jan 2023
classdef ExperimentMeta < dj.Part
    properties (SetAccess = protected)
        master = ns.Experiment;  % Part  table for the Experiment
    end

    methods
        function o = ExperimentMeta(varargin)
            o = o@dj.Part(varargin{:});
        end
    end

end