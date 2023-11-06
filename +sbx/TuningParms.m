%{
# Table with tuning analysis parameter sets.
tuningtag                   : varchar(128)                  # Short name to idenfity this tuning analysis
---
description=null            : varchar(1024)                 # Description of the parameter set
parms                       : longblob                      # struct of all parameters
%}
% 
% 
% 
% BK - Sept 2023
classdef TuningParms < dj.Lookup
    
end
