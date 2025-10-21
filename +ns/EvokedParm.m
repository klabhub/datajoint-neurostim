%{
# ns.EvokedParm lookup table for ns.Evoked
# group_fun isfunction to sort trials, if empty, it will default to group 
# by dimension

evtag : varchar(32) # unique tag
etag : varchar(32) # Epoch etags to select epochs to average across
---
group_fun = NULL : varchar(32)
epoch_query = 'flag = ""' : varchar(32)
%}

classdef EvokedParm < dj.Lookup & dj.DJInstance
    
end