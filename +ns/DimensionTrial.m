%{
# One row per trial in a condition - enables relational joins on trial number
-> ns.DimensionCondition
trial : int    # Trial number
%}
% This table expands the trials blob in ns.DimensionCondition into individual
% rows, allowing condition lookup via a relational join on trial number rather
% than blob membership tests (ismember).
%
% Example: get condition name for each EpochChannel row
%   result = fetch(ns.EpochChannel & restriction * ns.DimensionTrial, 'trial', 'channel', 'name');
%
% Example: restrict EpochChannel to a single condition
%   ec = ns.EpochChannel & (ns.DimensionTrial & 'name="gabor:orientation:30"');
%
% See also ns.DimensionCondition, ns.Dimension
% BK - May 2026.
classdef DimensionTrial < dj.Part
    properties (SetAccess = protected)
        master = ns.Dimension
    end
end
