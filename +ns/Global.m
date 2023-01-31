%{
# A table with global properties (e..g the data root)
id : smallint auto_increment
---
name: varchar(255) # Name
value : longblob  # Any value.
%}
%
% BK - April 2022

classdef Global < dj.Manual

    methods (Access=public)
        function v = get(tbl,prm)            
            tpl = tbl & struct('name',prm);
            if count(tpl)==1
                v= fetch1(tpl,'value');
            else
                v = '';
            end
        end
    end
end