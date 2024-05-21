%{
# A session refers to all recordings on a single day from a single subject
session_date: date  #Recording date (ISO 8601)
-> ns.Subject       #Subject (FK)
---
%}
% BK - April 2022
classdef Session < dj.Manual
    methods (Access=public)
        
        function v = folder(tbl)
             data = fetch(tbl,'session_date');
             dates = cellfun(@(x) datestr((x),'YYYY/mm/DD'),{data.session_date},'UniformOutput',false); %#ok<DATST>
             v = string(fullfile(getenv('NS_ROOT'),dates'));
        end
    end
end