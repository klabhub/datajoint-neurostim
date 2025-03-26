function [signal,time,channelInfo,recordingInfo] = read(key,parms)
% The key refers to an sbx file in an experiment.
%  This function will first check that preprocessed data exist for
% the parms.prep preprocessing instructions (Which should match a row in
% sbx.PreprocessedParm. An error is generated if no preprocessed data exist.
% If the preprocessed data exist, this funciton extracts (from files on disk)
% the relevant aspects to store in ns.C.
%
% The parms struct contains the following fields
% .prep -  A unique name to identify the preprocessing instructions (= a
% row in sbx.PreprocessedParm)
% .what - 'F','Fneu', or 'spks'
%
% To link  channels in ns.C with ROIs, run populate(sbx.Roi) after
% populating the ns.C table.
%

p = sbx.Preprocessed;
ks = p.keySource;
if ~exists(ks&key)
    fprintf('No files to analyze in this session\n');
    return;
end
if ~exists(sbx.Preprocessed & key & struct('prep', parms.prep))
    % The session has not been preprocessed. Do that first.
    error('No preprocessed data for %s in session %s for subject %s. Run populate(sbx.Preprocessed,prep="%s") first',parms.prep,key.session_date,key.subject,parms.prep);
end

% Read the npy results from the suite2p folder and store them in
% the table.
%% Determine nstime of each frame in this session
thisSession =(ns.Session & key);
allExptThisSession = ns.Experiment & (ns.File & 'extension=''.sbx''') &thisSession;
analyzeExptThisSession = analyze(allExptThisSession,strict=false);
nrFramesPrevious = 0;
for exptThisSession = fetch(analyzeExptThisSession,'ORDER BY starttime')'
    % Get the info structure that sbx saves
    info = sbx.readInfoFile(exptThisSession);
    nrFrames  = info.nrFrames;
    nrPlanes = info.nrPlanes;
    if strcmpi(exptThisSession.starttime,key.starttime)
        mdaq = proj(ns.C & 'ctag=''mdaq'''&exptThisSession,'time')* proj(ns.CChannel  & 'name=''laserOnDig''','signal');
        assert(exists(mdaq),'%s does not have the requred mdaq//laserOnDig channel yet. populate it first',exptThisSession.starttime)
        laserOnTTL = fetch(mdaq,'signal','time');
        laserOnIx = diff(laserOnTTL.signal)>0.5; % Transition from 0-1
        nstime = linspace(laserOnTTL.time(1),laserOnTTL.time(2),laserOnTTL.time(3));
        frameNsTime = nstime(laserOnIx);        % Time in ns time.

        % There always appears to be 1 extraneous TTL at the start
        frameNsTime(1) =[];
        nrTTL = numel(frameNsTime);


        % Sanity check
        delta = nrFrames- floor(nrTTL/nrPlanes) ;
        % Allow a slack of 3 ttls. TODO: make a prep parameter.
        if delta >0 || delta <=-3
            error('Cannot map SBX frames to trials; TTL-Frame mismatch (%d TTL %d frames in sbx).\n',nrTTL,nrFrames);
        else
            % Assume there were additional extraneous TTLs at the start.
            frameNsTime=frameNsTime(-delta+1:end);
        end

        keepFrameIx = nrFramesPrevious+(1:nrFrames);
        break; % We have what we need; break the loop over experiments in this session
    else
        nrFramesPrevious = nrFramesPrevious +nrFrames;
    end
end



%% Read the .npy or .mat Output
fldr= fullfile(folder(ns.Experiment & key),fetch1(sbx.Preprocessed & key & struct('prep',parms.prep),'folder'));
planes = dir(fullfile(fldr,'plane*'));
recordingInfo =sbx.readInfoFile(key);  % Store the info struct
signal=[];
rois = [];
maxRoi = 0;
assert(~isempty(planes),"No preprocessed files in %s ",fldr);
for pl = 1:numel(planes)
    %% Read npy
    tic;
    fprintf('Reading numpy files...\n')
    thisFile = fullfile(fldr,planes(pl).name,[parms.what '.npy']);
    if ~exist(thisFile,"file")
        error('File %s does not exist',thisFile);
    end
    thisSignal = single(py.numpy.load(thisFile,allow_pickle=true));
    thisSignal = thisSignal(:,keepFrameIx)';
    signal = [signal  thisSignal]; %#ok<AGROW>
    fprintf('Done in %s.\n',seconds(toc))

end
[nrFrames,nrROIs] = size(signal);
time = [frameNsTime(1) frameNsTime(end) nrFrames];
% Note that ROI are numbered across planes. This is matched in
% sbx.PreprocessedRoi to allow inner joins with CChannel.
% (see example in sxb.PreprocessedRoi/plotSpatial)
channelInfo =  struct('nr',num2cell(1:nrROIs)');
end
