%{ Skeleton Pair
%-> (source) -> ns.CChannel(channel)
%-> (target) -> ns.CChannel(channel)
%---
%
%}
% 

classdef SkeletonPair <dj.Part
     properties (SetAccess=protected)
        master = fc.Skeleton;
     end
end