%{
# A condition refers to a set of trials in an Experiment with matching parameters.
-> ns.Dimension
name        : varchar(128)              # Condition name (plg_prm_value)
---
trials : blob                            # The trials in this condition
%}
% See also ns.Dimension
% BK - March 2023.
classdef DimensionCondition < dj.Part
     properties (SetAccess = protected)
        master = ns.Dimension
    end
end