function [signal,time,channelInfo,recordingInfo] = read(key,parms)
% The key refers to an sbx file in an experiment.
%  This function will preprocess the raw data  (for the whole session),
% create a sbx.Preprocessed table and then
% extract the relevant aspects to store in ns.C. Preprocessing is not
% repeated if the Preprocessed data already exist. I.e., in a session
% with multiple experiments, preprocessing is done when ns.C is populated
% for the first experiment and subsequent calls to populate read the
% output files and add relevant content to ns.C.
%
% The parms struct contains the following fields
% .toolbox - 'suite2p', or 'caiman'
% .ops -  A struct with options passed to the toolboxes.
% .prep -  A unique name to identify these preprocessing instructions.
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
    make(sbx.Preprocessed,ns.stripToPrimary(sbx.Preprocessed,key),parms)
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



%% Read the NPY Output
fldr= fullfile(folder(ns.Experiment & key),fetch1(sbx.Preprocessed & key & struct('prep',parms.prep),'folder'));
planes = dir(fullfile(fldr,'plane*'));
recordingInfo =sbx.readInfoFile(key);  % Store the info struct
signal=[];
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
[nrFrames,nrROIs] = size(signal); %#ok<ASGLU>
time = [frameNsTime(1) frameNsTime(end) nrFrames];

channelInfo =  struct('nr',num2cell(1:nrROIs)');
 
end
