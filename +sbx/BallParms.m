%{
#  Parameter sets used for ball tracking 
    tag:  varchar(128)
    ---
    description : varchar(1024)
    parms : longblob  # struct of all parameters
%}
% 
% EXAMPLE
% The user fills this table with named parameter sets of their chosing that
% they then use to populate the Ball table. 
% 
% Create a struct with parameters for a toolbox and add this to the table. 
% For instance, a parameter setting to use for Ball tracking recordings and
% the xcorr algorithms implemented in sbx.Ball:
%
% xParms = struct;  
% add this to this table by calling: 
% insert(sbx.BallParms, struct('tag','xcorr', ...
%                                 'description', 'Using xcorr to estimate ball motion',...
%                                 'parms',xParms))
% 
% Then populate(sbx.Ball,'tag=''xcorr'') will populate the
% Ball table with ball tracking data based on these settings.
% 
% BK - Sept 2023
classdef BallParms < dj.Lookup
    
end