%{
#  Table with tuning analysis parameter sets.  
    tuningtag:  varchar(128)    # Short name to idenfity this tuning analysis 
    ---    
    description  = NULL : varchar(1024)  # Description of the parameter set 
    parms : longblob  # struct of all parameters
%}
% 
% For example, see sbx.Tuning 
% 
% BK - Sept 2023
classdef TuningParms < dj.Lookup
    
end