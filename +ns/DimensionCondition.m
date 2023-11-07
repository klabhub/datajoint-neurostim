%{
# A condition refers to a set of trials in an Experiment with matching parameters.
-> ns.Dimension
name        : varchar(128)      # Condition name (plg_prm_value)
---
value       : blob              # The value of the parameters for this condition
trials      : blob              # The trials in this condition
%}
% See also ns.Dimension
% BK - March 2023.
classdef DimensionCondition < dj.Part
     properties (SetAccess = protected)
        master = ns.Dimension
    end
end