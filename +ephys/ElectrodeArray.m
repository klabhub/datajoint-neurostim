%{
# A table of permanently implanted electrode arrays or EEG headcaps 
arrayname : varchar(255) #  A unique name identifying the array
---
nrelectrodes : int        # Number of electrodes in this array
description = NULL : varchar(1024) #  A brief description
%}
%

classdef ElectrodeArray< dj.Lookup
    
  
end