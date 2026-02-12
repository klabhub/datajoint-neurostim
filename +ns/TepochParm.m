%{
# TepochParm: parameters to transform epoch data.
ttag : varchar(32)  # tag for this transformed epoch
etag : varchar(32) # which epochs to transform 
---
fun : varchar(32) # Function to do the transform (see ns.cache/compute)
options: blob     # Options passed to the compute function
window  = NULL  : tinyblob          # Start and stop time of the epoch.Defaults to entire epoch
channels =NULL  : blob              # Channels to include. Defaults to all in the etag
grouping = NULL :blob               # Grouping to use. Defaults to no grouping.          
trial   = NULL : blob               # Trials to include. Defaults to all in the etag    
%}


classdef TepochParm < dj.Lookup & dj.DJInstance

end