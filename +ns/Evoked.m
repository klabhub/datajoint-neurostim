%{
# ns.Evoked is average ERPs per channel and dimension
# by default, grouped by dimension conditions
# if a group_fun provided in ns.EvokedParm, it is used to further split
# into groups.
# Part tables, ns.EvokedChannel contain actual signal
-> ns.EvokedParm
-> ns.DimensionCondition
group : varchar(32)
---
n_epoch: int
%}

classdef Evoked < dj.Computed & dj.DJInstance

    properties

        keySource
        trials_ dj.DJProperty = dj.DJProperty.empty

    end

    properties (Dependent)

        trials

    end
    

    methods

        function v = get.keySource(varargin)

            % only those events with some clean epochs
            v = ns.DimensionCondition * (ns.EvokedParm) & (ns.Epoch & 'flag=""');

        end


    end


    methods (Access = protected)

        function makeTuples(evTbl, key)

            evTbl


        end

    end

end