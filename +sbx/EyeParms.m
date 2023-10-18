%{
#  Parameter sets used for eye position and pupil size extraction
    tag:  varchar(128)
    ---
    description : varchar(1024)
    parms : longblob  # struct of all parameters
%}
% 
% EXAMPLE
% The user fills this table with named parameter sets of their chosing that
% they then use to populate the Eye table. 
% 
% Create a struct with parameters for a toolbox and add this to the table. 
% For instance, a parameter setting to use for eye tracking recordings and
% the imfindcircles algorithms implemented in ns.Pose:
% (see help imfindcircles for the meaning of these parameters).
% imfindCirclesParms = struct(''ObjectPolarity','bright',...
%                                   'Method','PhaseCode',...
%                                   'Sensitivity',0.85,...
%                                    'EdgeThresdhold',[],...
%                                     'radiusRange',[6 20])
%                     
% add this to this table by calling: 
% insert(sbx.EyeParms, struct('tag','imfindcircles', ...
%                                 'description', 'Using imfindcircles to detect pupil location and area', ...
%                                 'parms',imfindcirclesParms))
% 
% Then populate(sbx.Eye,'tag=''imfindcircles'') will populate the
% Pose table with eye tracking data based on these settings.
% 
% BK - Sept 2023
classdef EyeParms < dj.Lookup
    
end