%{
# A table of permanently implanted electrode arrays or EEG headcaps 
array : varchar(255) #  A unique name identifying the array
---
channels : blob           # The channesl recording from the electrodes in this array
description = NULL : varchar(1024) #  A brief description
%}
%

classdef Array< dj.Lookup
    
  
end