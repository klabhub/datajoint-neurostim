% -------------------------------------------------------------------------
% ns.SessionInfoParms
% -------------------------------------------------------------------------
%{
# ns.SessionInfoParms: Parameters for defining session ordering
sinf_tag : varchar(32)
ctag : varchar(32) #primary
paradigm : varchar(64)
---
%}
classdef SessionInfoParms < dj.Lookup & dj.DJInstance
end