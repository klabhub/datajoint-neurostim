%{
# Functional Connectivity between channels
-> fc.Fc
source    : int  # Channel number of the source of the FC
target    : int  # Channel number of the target of the FC
---
fc : float # The functional connectivity between the source and the target
p  = NULL : float # p-value associated with the FC.
err = NULL : float # a measure of uncertainty associated with the FC.
%}

classdef FcPair < dj.Part
     properties (SetAccess=protected)
        master = fc.Fc;
     end
    % Insertion is done by the master fc.Fc/makeTuples
end