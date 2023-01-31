%{
# Meta data for a Neurostim subject. Read from a JSON file.
-> ns.Subject
meta_name : varchar(255)     # The name of the meta data
---
meta_value : varchar(255)   # The value of the meta data
%}
% BK - Jan 2023
classdef SubjectMeta < dj.Part
    properties (SetAccess=protected)
        master =ns.Subject;
    end
    
end
