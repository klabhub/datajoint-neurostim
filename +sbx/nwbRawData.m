function nwbRawData(fileTpl,nwbRoot,pv)
% Function that reads an sbx data file, and adds it to the nwbTrgt (an
% instance of NwbFile). This is meant to be called automatically from
% nwbExport.
% Values that cannot be read from the info file are read from
% pv.passthrough.sbx (or given a hardcoded default value here)
%
%
arguments
    fileTpl (1,1) struct
    nwbRoot (1,1) NwbFile
    pv (1,1) struct
end

if ~isfield(pv.passthrough,'sbx')
    error('Please provide sbx information in the call to nwbExport')
end

info = sbx.readInfoFile(ns.Experiment & fileTpl);
device = types.core.Device('manufacturer','NeuroLabWare','description','TwoPhoton Microscope');
nwbRoot.general_devices.set('Device',device);
optical_channel = types.core.OpticalChannel( ...
    'description', 'description', ...
    'emission_lambda',getParms(pv.passthrough.sbx,'emission',500));
imaging_plane = types.core.ImagingPlane('optical_channel',optical_channel, ...
    'device',types.untyped.SoftLink(device), ...
    'excitation_lambda',getParms(pv.passthrough.sbx,'excitation',510), ...
    'imaging_rate',getParms(pv.passthrough.sbx,'framerate',15.5), ...
    'indicator',getParms(pv.passthrough.sbx,'indicator','GCamp6s'),...
    'location',getParms(pv.passthrough.sbx,'location','V1'), ...
    'description',getParms(pv.passthrough.sbx,'description','SBX Imaging Plane'));

nwbRoot.general_optophysiology.set('imaging_plane',imaging_plane);


ff = strrep(fullfile(folder(ns.Experiment & fileTpl),fileTpl.filename),'.sbx','');

if getParms(pv.passthrough.sbx,'meanOnly',true)
    % Instead of the raw data, store only the mean image. 
    data = types.untyped.DataPipe(permute(fetch1(sbx.Preprocessed & fileTpl,'img'),[3 1 2])); % Compress for consistencey with raw image data.
else
    % Read the raw image data. This works, but the compression (an nwb save
    % time) seems to take forever nd for a 30GB file needs more thatn 100GB
    % of RAM? TODO: troubleshooting
    fprintf('Reading %s ...',fileTpl.filename); tic
    data = sbxread(char(ff),0,info.nrFrames,info);
    fprintf('Done (%s s).\n',seconds(toc))
end
% TODO: the timing of the raw images can be extracted om the ns.C
% that stores a ROI time series. 
tpi = types.core.TwoPhotonSeries( ...
    'imaging_plane',types.untyped.SoftLink(imaging_plane), ...
    'starting_time',0.0, ...
    'starting_time_rate',getParms(pv.passthrough.sbx,'framerate',15.5), ...
    'starting_time_unit','seconds',...
    'data',data,...
    'data_unit','lumens', ...
    'data_continuity','continuous', ...
    'description','Two photon data',...
    'field_of_view',([info.dxcal info.dycal].*info.sz)*1e-6,...
    'format','raw',...
    'pmt_gain',info.config.pmt0_gain...
    );
nwbRoot.acquisition.set('2pInternal',tpi);
end

function v = getParms(s,fld,default)
% Try to read s.fld, if it does not exist, warn and return default
if isfield(s,fld)
    v =s.(fld);
else
    fprintf('sbx.%s does not exist. Using default value: %s\n',fld,string(default));
    v = default;
end
end



function x = sbxread(fname,k,N,info)
% img = sbxread(fname,k,N,varargin)
%
% Reads from frame k to k+N-1 in file fname
%
% fname - the file name (e.g., 'xx0_000_001')
% k     - the index of the first frame to be read.  The first index is 0.
% N     - the number of consecutive frames to read starting with k.
%
% If N>1 it returns a 4D array of size = [N rows cols #pmt]
% If N=1 it returns a 3D array of size = [rows cols #pmt]
%
% #pmts is the number of pmt channels being sampled (1 or 2)
% rows is the number of lines in the image
% cols is the number of pixels in each line
% N is the number of frames
%
% Adapted from readsbx by Dario Ringach


if(info.scanmode==0)
    info.recordsPerBuffer = info.recordsPerBuffer*2;
end

if(~isfield(info,'chan'))
    switch info.channels
        case 1
            info.nchan = 2;      % both PMT0 & 1
            factor = 1;
        case 2
            info.nchan = 1;      % PMT 0
            factor = 2;
        case 3
            info.nchan = 1;      % PMT 1
            factor = 2;
    end
else
    info.nchan = info.chan.nchan;
end

% Open the sbx file
fid = fopen([fname '.sbx']);
d = dir([fname '.sbx']);
info.nsamples = (info.sz(2) * info.recordsPerBuffer * 2 * info.nchan);   % bytes per record
if isfield(info,'scanbox_version')
    switch(info.scanbox_version)
        case 2
            info.max_idx =  d.bytes/info.recordsPerBuffer/info.sz(2)*factor/4 - 1;
            info.nsamples = (info.sz(2) * info.recordsPerBuffer * 2 * info.nchan);   % bytes per record
        case 3
            info.max_idx = d.bytes/prod(info.sz)/info.nchan/2 -1;
            info.nsamples = prod(info.sz)*info.nchan*2;
        otherwise
            error('Invalid Scanbox version');
    end
else
    info.max_idx =  d.bytes/info.bytesPerBuffer*factor - 1;
end

% Read all data.
try
    fseek(fid,k*info.nsamples,'bof');
    x = fread(fid,info.nsamples/2 * N,'uint16=>uint16');
    x = reshape(x,[info.nchan info.sz(2) info.recordsPerBuffer  N]);
catch
    error('Cannot read frame.  Index range likely outside of bounds.');
end
% Reorder to match NBW requirement of time first.
x = intmax('uint16')-permute(x,[4 3 2 1]); % TODO: frame first or last >?Change order to have [Frames Cols Rows Channels]

% folding lines

if isfield(info,'fold_lines')
    if info.fold_lines>0
        error('fold lines - not implemented yet')
        % sx = size(x);
        % x = reshape(x,sx(1),info.fold_lines,[],sx(3),sx(4));
        % x = permute(x,[1 2 4 3 5]);
        % sx = size(x);
        % x = reshape(x,sx(1),sx(2),sx(3),[]);

        % %             t = info.frame(idx)*512+info.line(idx);   % recode TTLs
        % %             info.frame = floor(t/info.fold_lines);
        % %             info.line = mod(t,info.fold_lines);

    end

end

fclose(fid);

% Setup for compression. 
x = types.untyped.DataPipe('data',x); % Not clear how to choose ,'chunkSize', [info.sz));

end
