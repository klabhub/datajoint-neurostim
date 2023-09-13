%{
#  Preprocessing parameter sets used for calcium imaging data
    tag:  varchar(128)
    ---
    toolbox : enum('suite2p','caiman')
    description : varchar(1024)
    parms : longblob  # struct of all parameters
%}
% 
% EXAMPLE
% The user fills this table with named parameter sets of their chosing that
% they then use to populate the Preprocessed table. 
% 
% Create a struct with parameters for a toolbox and add this to the table. 
% For instance, a parameter setting to use for Scanbox recordings with gGamp6s :
% gcamp6s = struct('fs',15.5, ...    % Framerate on ScanBox/Neurolabware TPI
%                      'tau',1.3,... % Recommended for Gcamp6s
%                      'delete_bin',true, ...   % No need to keep
%                      'use_builtin_classifier',true); % Use the built-in classifier
% add this to this table by calling: 
% insert(sbx.PreParameters, struct('tag','gCamp6s', ...
%                                 'toolbox','suite2p', ...
%                                 'description', 'Longer tau for GCamp6s, deleting bin file, builtin classifier', ...
%                                 'parms',gCamp6s))
% 
% Then populate(sbx.Preprocessd,'prep=''gcamp6s'') will populate the
% Preprocessed table with all default suite2p settings, except those set in
% gGamp6s struct above.
% 
% BK - Mar 2023
classdef PrepParms < dj.Lookup
    
end