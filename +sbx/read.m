function [signal,time,channelInfo,recordingInfo] = read(key,parms)
% The key refers to an sbx file in an experiment.
%
%  This function will first check that preprocessed data exist for
% the parms.prep preprocessing instructions (Which should match a row in
% sbx.PreprocessedParm. An error is generated if no preprocessed data exist.
%
% Once the preprocessed data exist, this function extracts (from files on disk)
% the relevant aspects to store in ns.C.
%
% The parms struct contains the following fields
% .prep -  A unique name to identify the preprocessing instructions (= a row in sbx.PreprocessedParm)
% .what - 'F'    - Fluorescence
%         'Fneu'  - Neuropil
%         'spks'  - deconvolved spikes from suite2p's OASIS
%         
% .neuropilfactor - If .what ="F" then setting this to a non-zero value will compute
%               bacground corrected fluorescence (F-factor*Fneu)
% .restrict - A restriction on sbx.PreprocessedRoi to limit the ROI to a subset.
%               For instance 'pcell>0.75'

prep = sbx.Preprocessed;
ks = prep.keySource;
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
[keepFrameIx,frameNsTime] = sbx.framesForExperiment(key);
      

%% Read the .npy or .mat Output
fldr= fullfile(folder(ns.Experiment & key),fetch1(sbx.Preprocessed & key & struct('prep',parms.prep),'folder'));
planes = dir(fullfile(fldr,'plane*'));

signal=[];
rois = [];
maxRoi = 0;
for pl = 1:numel(planes)
    %% Read npy
    tic;
    if ismember(parms.what, ["F" "Fneu" "spks"])
        % Suite 2p numpy files
        fprintf('Reading numpy files...\n')
        thisFile = fullfile(fldr,planes(pl).name,[parms.what '.npy']);
        if ~exist(thisFile,"file")
            error('File %s does not exist',thisFile);
        end
        thisSignal =  ndarrayToArray(py.numpy.load(thisFile,allow_pickle=true),single=true);
        framesNotInFile = sum(keepFrameIx > size(thisSignal,2));
        assert(framesNotInFile==0,"%s has %d too few frames for %s on %s",thisFile,framesNotInFile,key.starttime, key.session_date);
        thisSignal = thisSignal(:,keepFrameIx)';
        if  parms.what=="F" && isfield(parms,'neuropilfactor') && parms.neuropilfactor ~=0
            % Compute background/neuropil corrected fluorescence
            neuFile = strrep(thisFile,'F.npy','Fneu.npy');
            Fneu  =  ndarrayToArray(py.numpy.load(neuFile,allow_pickle=true),single=true);
            Fneu = Fneu(:,keepFrameIx)';
            thisSignal = thisSignal - parms.neuropilfactor*Fneu;
        end
        rois = [rois ;(1:size(thisSignal,2))'+maxRoi]; %#ok<AGROW>
        maxRoi = max(rois);
    else
        error('Unknown what (%s) for sbx.read ',parms.what);
    end
    signal = [signal  thisSignal]; %#ok<AGROW>
    fprintf('Done in %s.\n',seconds(toc))

end

if isfield(parms,'restrict')
    % Restrict with a query on sbx.PreprocessedRoi
    keepRoi = [fetch((sbx.PreprocessedRoi & parms.restrict ) & key ,'roi').roi];
    signal = signal(:,keepRoi);
    rois    = keepRoi;
end
nrFrames = size(signal,1);
time = [frameNsTime(1) frameNsTime(end) nrFrames];
% Note that ROI are numbered across planes. This is matched in
% sbx.PreprocessedRoi to allow inner joins with CChannel.
% (see example in sxb.PreprocessedRoi/plotSpatial)
channelInfo =  struct('nr',num2cell(rois)');
recordingInfo =struct('dummy',true); 
end