%{ 
# Functional Connectivity skeleton 
-> ns.Session
-> ns.FcSkeletonParm  # Parameters determining how to estimate the skeleton
---
source : blob                 #
target: blob
%}

classdef FcSkeleton <dj.Computed

end