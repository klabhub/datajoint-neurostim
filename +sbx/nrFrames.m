function v = nrFrames(sbxFile,info)
% Determine the number of frames in an sbx file.
% INPUT
% sbxFile = full path to an sbxFile
% OUTPUT
% N = number of frames in the file.
% 
arguments
    sbxFile (1,1) string
    info (1,1) struct
end

if ~exist(sbxFile,"file")
    error('This file does not exist: %s',sbxFile)
end

[fldr,file,~] =fileparts(sbxFile);
if isempty(info)
    load(fullfile(fldr,file +".mat"),'info');
end
dirInfo = dir(fullfile(fldr,file + ".sbx"));
% Determine the number off planes  from the optotune settings.
if isfield(info,'otwave')
    nrPlanes = max(1,numel(info.otwave));
elseif isnan(info.otparam(3))
    nrPlanes = 1 ;
else
    nrPlanes = info.otparam(3);
end
nrChannels = numel(info.channels);
% 2 bytes per pixel.
v= dirInfo.bytes./prod(info.sz)/nrChannels/nrPlanes/2;