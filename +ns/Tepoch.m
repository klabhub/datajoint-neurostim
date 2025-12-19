%{
# Transformed Epoch. 
-> ns.Epoch
-> ns.TepochParm
---
x : longblob
description 
%}
classdef Tepoch < dj.Computed & dj.DJInstance   
    methods (Access = protected)

        function makeTuples(tbl, key)

            
        end
    end

end