% -------------------------------------------------------------------------
% ns.EpochInfoData
% -------------------------------------------------------------------------
%{
# ns.EpochInfoData: True order of valid epochs
-> ns.EpochInfo
-> ns.Epoch
---
trial_order=null : int  # Sequential order, NULL if excluded
%}
classdef EpochInfoData < dj.Part
    properties(SetAccess=protected)
        master = ns.EpochInfo
    end
end