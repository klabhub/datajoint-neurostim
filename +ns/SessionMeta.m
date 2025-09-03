%{
# Meta data for a Neurostim Session. Read from a JSON file.
-> ns.Session
meta_name : varchar(255)     # The name of the meta data
---
meta_value : varchar(2048)  # The value of the meta data
%}
% nsScan reads from JSON and nsAddToDataJoint inserts this information into
% this table. Use the nsMeta app to update meta data.
% BK - Jan 2023
classdef SessionMeta < dj.Part
    properties (SetAccess=protected)
        master =ns.Session;
    end
end
