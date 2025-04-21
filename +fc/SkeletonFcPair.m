%{
# Skeleton FC 
-> fc.Skeleton
source    : int  # Channel number of the source of the skeleton FC
target    : int  # Channel number of the target of the skeleton FC
---
fc : float # The functional connectivity between the source and the target
p  = NULL : float # p-value associated with the FC.
err = NULL : float # a measure of uncertainty associated with the FC.
%}

classdef SkeletonFcPair < dj.Part
     properties (SetAccess=protected)
        master = fc.Skeleton;
     end
    % Insertion is done by the master fc.Skeleton/makeTuples
end