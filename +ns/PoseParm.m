%{
#  Parameter sets used for DLC pose extraction
    tag:  varchar(128)
    ---
    description : varchar(1024)
    parms : longblob  # struct of all parameters used in the DLC process
%}
% 
% EXAMPLE
%    .config
%    .shuffle (1,1) double {mustBeNonnegative,mustBeInteger}
%    .trainingsetindex (1,1) double {mustBeNonnegative,mustBeInteger}
%    .TFGPUinference  (1,1) logical
         
% BK - Sept 2023
classdef PoseParm < dj.Lookup
    
end