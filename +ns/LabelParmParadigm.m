%{
# List of paradigms to which a LabelParm definition applies.
->ns.LabelParm
->ns.Paradigm 
---
%}
% Note that this table is used to restrict the application of LabelParm entries to specific paradigms.
% If a LabelParm entry has no corresponding rows in LabelParmParadigm, it applies to all paradigms.
classdef LabelParmParadigm < dj.Part
    properties (SetAccess = protected)
        master = ns.LabelParm
    end
end