%{
# Functional Connectivity between channels
(source)-> ns.CChannel(channel)  # The source channel
(target)-> ns.CChannel(channel)  # The target channel
-> ns.Fc
---
fc : float # The functional connectivity between the source and the target
p  = NULL : float # p-value associated with the FC.
err =NULL : float # a measure of uncertainty associated with the FC.
%}

classdef Pair < dj.Part
     properties (SetAccess=protected)
        master = ns.Fc;
     end
    % Insertion is done by the master ns.Fc/makeTuples
end