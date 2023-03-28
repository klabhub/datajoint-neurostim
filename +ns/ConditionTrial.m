%{
# Trials per ns.Condition
->ns.Condition
trial :smallint unsigned
%}
% See ns.Condition
classdef ConditionTrial < dj.Part  
    properties (SetAccess = protected)
        master = ns.Condition
    end
    % This table is filled by the populate function of the 
    % ns.Condition master class.
end