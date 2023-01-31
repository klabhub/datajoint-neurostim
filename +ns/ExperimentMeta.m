%{
# Meta data for a Neurostim experiment. Read from a JSON file.
-> ns.Experiment
meta_name : varchar(255)     # The name of the meta data
---
meta_value : varchar(255)   # The value of the meta data
%}
% BK - Jan 2023
classdef ExperimentMeta < dj.Part
    properties (SetAccess = protected)
        master = ns.Experiment;  % Part  table for the Experiment
    end
    methods (Access= public)
        function updateFromJson(tbl,experimentKey)     
            if ~isstruct(json)
                if isempty(json)
                    json= regexprep(fetch1(ns.Experiment & experimentKey,'file'),['(\.mat$)'],'.json'); % Swap extension
                end
                if exist(json,"file")
                    json  = readJson(json);
                else
                    %No json file, nothing to do.
                    return;
                end
            end
            values= struct2cell(json);
            for i=1:numel(values)
                if ~(ischar(values{i}) || isstring(values{i}))
                    values{i} =jsonencode(values{i},"PrettyPrint",true);
                end
            end
            newTuples = struct('starttime',experimentKey.starttime, ...
                            'session_date',experimentKey.session_date, ...
                            'subject',experimentKey.subject,...
                            'meta_name', fieldnames(json),...
                            'meta_value', values);
           
            delQuick(tbl & experimentKey); % Delete old. 
            tbl.insert(newTuples);         % Put new.   
        end
    end
end