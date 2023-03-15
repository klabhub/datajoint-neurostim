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
nrPlanes = max(1,numel(info.otwave));
nrChannels = numel(info.channels);
% 2 bytes per pixel.
v= dirInfo.bytes./prod(info.sz)/nrChannels/nrPlanes/2;