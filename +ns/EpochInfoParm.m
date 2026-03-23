% -------------------------------------------------------------------------
% ns.EpochInfoParm
% -------------------------------------------------------------------------
%{
# ns.EpochInfoParm: Parameters for determining epoch order and exclusion
ei_tag : varchar(32)
etag : varchar(32)
---
exclusion_query='' : varchar(255) # SQL restriction string for ns.Epoch (e.g. 'flag!="badByFixation"')
%}
classdef EpochInfoParm < dj.Lookup
end