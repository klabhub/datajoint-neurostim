%{
# TepochParm: parameters to transform epoch data.
ttag : varchar(32)
etag : varchar(32) # which epochs to transform 
---
fun : varchar(32) # Function to do the transform (see ns.EpochChannel/compute)
options: blob     # Options passed to the compute function
%}


classdef TepochParm < dj.Lookup & dj.DJInstance

end