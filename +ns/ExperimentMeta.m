%{
# Meta data for a Neurostim experiment. 
-> ns.Experiment
meta_name : varchar(255)     # The name of the meta data
---
meta_value : longblob       # The value of the meta data
%}
% nsScan reads from JSON and nsAddToDataJoint inserts this information into
% this table. Use the nsMeta app to update meta data.
% 
% Unlike the meta data for Subjects and Sessions, here the values are blobs to store
% more complex (non-char) meta data. See sbx.addExperimentMeta for an example.
% BK - Jan 2023
classdef ExperimentMeta < dj.Part
    properties (SetAccess = protected)
        master = ns.Experiment;  % Part  table for the Experiment
    end   
end