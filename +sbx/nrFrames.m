function v = nrFrames(expt,info)
% Determine the number of frames in an sbx experiment .
% INPUT
% expt = Tuple from the ns.Experiment table (corresponding to an sbx experiment)
% OUTPUT
% N = number of frames in the file.
% 
arguments
    expt (1,1) struct
    info struct = struct([])
end

sbxFile = ((ns.File & expt) & 'extension=''.sbx''');
if count(sbxFile) ~=1
    error('nrFrames requires a single matching experiment tpl as input, not %d ',count(sbxFile))
end
fname = fetch1(sbxFile & expt,'filename');
filename  = folder(ns.Experiment & expt) + fname;
dirInfo = dir(filename);
if isempty(info)
    % Read the corresponding info struct
    [fldr,file,~] =fileparts(strrep(filename,'.sbx','.mat')); 
    load(fullfile(fldr,file +".mat"),'info');
end

% Determine the number of planes  from the optotune settings.
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