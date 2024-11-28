%{
#  Preprocessing instructions for SBX files
prep         :  varchar(255)     # A  unique name for these preprocessing instructions
---
parms       : longblob          # struct containing all parameters used by packages such as suite2p.
%}
%
% EXAMPLE
% For suite2p, the parms should be a struct with field toolbox =suite2p,
% and a field ops that defines suite2p options:
% struct('prep','gcamp6s','parms',
%               struct('toolbox','suite2p',
%                       'ops',struct(...
%                       'fs',15.5, ...    % Framerate on ScanBox/Neurolabware TPI
%                      'tau',1.3,... % Recommended for Gcamp6s
%                       'delete_bin',true, ...   
%                       'use_builtin_classifier',true,... % Use the built-in classifier
%                        'look_one_level_down',true)));% Scanbox 3 makes subsub folders
%
% BK - June 2023

classdef PreprocessedParm < dj.Lookup

end