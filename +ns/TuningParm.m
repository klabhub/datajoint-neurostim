%{
# Table with tuning analysis parameter sets.
-> ns.CParm
tuningtag                   : varchar(128)                  # Short name to idenfity this tuning analysis
---
description=null            : varchar(1024)                 # Description of the parameter set
dimension                    :varchar(255)                  # The name of the dimension for which this tuning is to be computed (see ns.Dimension)
parms                       : longblob                      # struct of all parameters used for the tuning curve estimation
%}
% 
% BK - Sept 2023
classdef TuningParm < dj.Lookup
    
end
