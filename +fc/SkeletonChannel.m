%{ 
# Skeleton Channels
-> fc.Skeleton
channel : int 
---
%}
 

classdef SkeletonChannel <dj.Part
     properties (SetAccess=protected)
        master = fc.Skeleton;
     end
end