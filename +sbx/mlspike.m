function [signal,time,channelInfo,recordingInfo] = mlspike(key,parms)
% Function that uses the spikes and xplor repositories from Deneux et al.
% 2016 to run spike deconvolution.
%
% This function is called when populating ns.C with the appropriate
% parameters.
%
% With autocalibration defined it can be restricted to a subset of the samples
% to reduce the time spent on calibration. The sigma is always estimated per roi. 
% If autocalibration is requested and fails, the ROI will not be added to
% the C table. 
%
% EXAMPLE
% Setup the parameters for the alogrithm (see spikes/spk_demo.m)
%
%parms.mlparms = sbx.mlspikeDefaults("DENEUX16","gcamp6s");
%parms.mlparms.algo.nspikemax =4;  % Allow 4 spikes per bin (15 Hz)
%parms.mlparms.dographsummary = false;
%parms.mlparms.display = 'none';
%
%parms.autocal.amin = 0.05;
%parms.autocal.amax =0.2;
%parms.autocal.taumin = 0.25;
%parms.autocal.taumax = 2;
%parms.autocal.maxamp =4  % A maxamp of 4 seems necessary in our data.
%parms.autocal.display = 'none'; % No figures
%
% parms.secsForCal = 300; % Use 300 seconds at the beginning, middle, and end for autocal.
% parms.prep = 'gcamp6s'
% parms.restrict = 'pcell>0.75 AND radius>2.185'
% parms.neuropilFactor = 0.7; % Substract neuropil with this factor
%
% Create a ns.CParm struct that uses these parms
% deneux = struct('ctag','deneux16',...       % Name of this C
%                   'description','Deneux spikes',...
%                    'extension','.sbx',...
%                    'fun','sbx.mlspike',...
%                    'parms',parms);
%
% insertIfNew(ns.CParm,deneux);
% Start filling the table (restricted to sessions with preprocessed data)
% populate(ns.C,'ctag="deneux16"',sbx.Preprocessed)
arguments
    key (1,1) struct
    parms (1,:) struct
end

assert(exist('xplor.m','file'),"The xplor repository must be on the path for sbx.mlspike");
assert(exist('spk_est.m','file'),"The spikes repository must be on the path for sbx.mlspike");
prep = fetch(sbx.Preprocessed & key & struct('prep', parms.prep),'*');
assert(~isempty(prep),'No preprocessed data for %s in session %s for subject %s. Run populate(sbx.Preprocessed,prep="%s") first',parms.prep,key.session_date,key.subject,parms.prep);
warning('off','backtrace');
assert(isfield(parms,"mlparms"),'mlspike parameters must define mlparms.');

parms.mlparms.dt = 1./(prep.framerate/prep.nrplanes); % Match dt to framerate
if ~isfield(parms,'neuropilFactor')
    parms.neuropilFactor = 0.7; % Default neuropil correction factor
end
if ~isfield(parms,'secsForCal')
    parms.secsForCal  = inf; % Calibrate on the entire timecourse
end
if ~isfield(parms,'allowNegative')
    parms.allowNegative = 0.05; % 5% negative samples allowed by default
end
%% Session based
% Check what has already been done
fldr= fullfile(folder(ns.Experiment & key),fetch1(sbx.Preprocessed & key & struct('prep',parms.prep),'folder'));
roiToDo = table;
for pl=0:prep.nrplanes-1
    if isfield(parms,'restrict')
        % Restrict with a query on sbx.PreprocessedRoi
        roi = [fetch((sbx.PreprocessedRoi & parms.restrict & struct('plane',pl)) & key ,'roi').roi]';
    end
    nrRoiThisPlane  = numel(roi);
    % Results will be saved in this file:
    filename = fullfile(fldr,"plane" + string(pl), string(roi) + "." + key.ctag + ".mlspike.mat");
    done = false(nrRoiThisPlane,1);
    for i=1:nrRoiThisPlane
        done(i) = exist(filename(i),"file");
    end
    plane = repmat(pl,nrRoiThisPlane,1);
    roiToDo = [roiToDo ; table(filename,done,roi,plane)]; %#ok<AGROW>
end
if any(~roiToDo.done)
    % Run deconvolution for the ROI that do not have a file on disk yet.
    deconResults = loop(roiToDo(~roiToDo.done,:),parms,fldr);    
end

%% Check that all are now done
for i=1:height(roiToDo)
    roiToDo.done(i) = exist(roiToDo.filename(i),"file");
end
if any(~roiToDo.done)
    fprintf('Deconvolution failed on the following ROIs :\n')
    roiToDo(~roiToDo.done,:)
end
failFraction = mean(~roiToDo.done);
roiToDo = roiToDo(roiToDo.done,:); % Crop
%% Experiment specific
% Now files with spiketimes exist on disk. Load and extract the relevant frames for the
% current experiment to return to the ns.C caller.

[keepFrameIx,frameNsTime] = sbx.framesForExperiment(key);
nrFramesThisExpt = numel(keepFrameIx);
nrFramesInSession = prep.nrframesinsession; % Total frames across all experiments in the session
signal = nan(nrFramesThisExpt,height(roiToDo));
channelInfo =  struct('nr',num2cell(roiToDo.roi),'quality',NaN,'tau',NaN,'sigma',NaN,'neg',NaN,'nan',NaN,'autocal',NaN);
% Note that ROI are numbered across planes. This is matched in
% sbx.PreprocessedRoi to allow inner joins with CChannel.
% (see example in sxb.PreprocessedRoi/plotSpatial)
tic;
fprintf('Reading spikeml mat files...\n')
%% Read saved results    
for i = 1:height(roiToDo)
    load(roiToDo.filename(i),'spikeTimes','result');
    if ~isfield(result,'autocal')
        result.autocal = NaN;
    end
    channelInfo(i).quality = result.quality;
    channelInfo(i).tau= result.tau;
    channelInfo(i).a= result.a;
    channelInfo(i).neg= result.neg;
    channelInfo(i).nan = result.nan;
    channelInfo(i).autocal = result.autocal;
    
    % Convert spike times to a signal with counts
    % spikeTimes are for the ENTIRE SESSION, so create time vector for all session frames
    thisSignal = brick.timevector(spikeTimes,(0:nrFramesInSession-1)*parms.mlparms.dt,'count');
    framesNotInFile = sum(keepFrameIx > numel(thisSignal));
    assert(framesNotInFile==0,"%s has %d too few frames for %s on %s",roiToDo.filename(i),framesNotInFile,key.starttime, key.session_date);
    % Extract only the frames for this specific experiment
    thisSignal = thisSignal(keepFrameIx);
    signal(:,i) = thisSignal;    
end
fprintf('Done in %s.\n',seconds(toc))
time = [frameNsTime(1) frameNsTime(end) nrFramesThisExpt];
% Store the mean quality and fail fraction for this experiment (session
% really)
mQuality = mean([channelInfo.quality],"omitmissing");
recordingInfo = struct('quality',mQuality,'fail',failFraction);
end

function [out] = loop(roiT,parms,fldr)
% Function to read the F/Fneu and start a pool of workers to run spikeml on
% each roi. Update messages are provided after each. 
arguments
    roiT table % List of rois that have not been done yet
    parms (1,1) struct
    fldr (1,1) string
end
global BRICKPROGRESS %#ok<GVMIS>
BRICKPROGRESS = false;

F = [];
for pl = unique(roiT.plane)
    thisFile = fullfile(fldr,"plane" +  string(pl) ,'F.npy');
    if ~exist(thisFile,"file")
        error('File %s does not exist',thisFile);
    end
    thisF=  ndarrayToArray(py.numpy.load(thisFile,allow_pickle=true),single=true);
    keepRoi = roiT.roi(roiT.plane==pl);

    thisFile = fullfile(fldr,"plane" +  string(pl) ,'Fneu.npy');
    if ~exist(thisFile,"file")
        error('File %s does not exist',thisFile);
    end
    thisFNeu=  ndarrayToArray(py.numpy.load(thisFile,allow_pickle=true),single=true);
    F = [F;(thisF(keepRoi,:)-parms.neuropilFactor*thisFNeu(keepRoi,:))']; %#ok<AGROW>
end
%% Use parallel pool if requested
pool = nsParPool;
[nrSamples, nrRoi] = size(F);
fprintf('Queue %d channels with %d samples for deconvolution at %s\n',nrRoi,nrSamples,datetime("now"))
tStart = tic;
nrRoiDone = 0;
for i=1:nrRoi
    future(i) = parfeval(pool,@deconvolve,1,F(:,i),parms,roiT.filename(i)); %#ok<AGROW>
end
afterEach(future,@afterDone,0,PassFuture = true);  % Update the command line
wait(future); % Wait for all to complete
keep = cellfun(@isempty,{future.Error});% Remove errors
assert(any(keep),'All %d ROI mlspike deconvolve failed in %s ',nrRoi,fldr);
out = fetchOutputs(future(keep));

warning('on','backtrace');

    function afterDone(ftr)
        % Nested message to command line that a job is done or errored out
        nrRoiDone = nrRoiDone +1;
        if ~isempty(ftr.Error)
            fprintf("%s failed after %s - (Error: %s)\n ",func2str(ftr.Function),ftr.RunningDuration,ftr.Error.message);
            if ~isempty(ftr.Diary)
                ftr.Diary
            end
        else
            cumTime= seconds(toc(tStart));
            cumTime.Format = "hh:mm:ss";
            secsPer = cumTime/nrRoiDone;
            timeRemaining = (nrRoi-nrRoiDone)*secsPer;
            fprintf("%s is done (%d/%d) after %s - (State: %s). Cumulative time: %s. Remaining: %s \n ",func2str(ftr.Function),nrRoiDone,nrRoi,ftr.RunningDuration,ftr.State,cumTime,timeRemaining);
        end
    end

end

function result = deconvolve(F,parms,filename)
% Function that does the actual mlspike deconvolution and saves the results
% to disk in a file named after the ctag. 
arguments
    F (:,1) single  % Fluorescence signal
    parms (1,1) struct   % parms
    filename  (1,1) string    
end

isNaN = isnan(F);
F(isNaN) = 0; % Could remove samples instead or linearly interpolate, but this should be rare (missing F)
isNegative  = F<0;
% F is supposed to be raw fluoresence and therefore positive, but the neuropil
% correction can occasionally generate negative F.  Set those to zero.
F(isNegative) = 0;

doAutoCal = isfield(parms,'autocal') ;

if mean(isNegative)>parms.allowNegative
    % Fail if more than 5% are negative
    result = struct('quality',0,'tau',nan,'a',nan,'sigma',nan,'neg',mean(isNegative),'nan',mean(isNaN),'autocal',nan);
    return;
else
    if any(isNegative)
        fprintf('%.2f %% of samples have negative F. Set to 0.\n',100*mean(isNegative));        
    end    
    if doAutoCal
        % Autocalibration with the specified parameters
        pax = spk_autocalibration('par'); % Get defaults, then overrule with parms.autocal
        % Copy values from parms.autocal struct to pax
        pax = brick.structmerge(pax,parms.autocal,'strict','recursive','type');
        pax.dt = parms.mlparms.dt;
        pax.mlspikepar = parms.mlparms;
        pax.mlspikepar.dographsummary = false;
        if isfinite(parms.secsForCal)
            % Use a subset of the F to determine calibration.
            % secsForCal at the start, middle, and end of the session
            nrSamples = size(F,1);
            t = (0:nrSamples-1)*pax.dt;
            useForCal = t < parms.secsForCal | t>(t(nrSamples)-parms.secsForCal) |  abs(t-t(round(nrSamples/2))) < 0.5*parms.secsForCal;
            [tauEst,ampEst,sigmaEst,events] = spk_autocalibration(F(useForCal),pax);
        else
            % Calibrate using the full F
            [tauEst,ampEst,sigmaEst,events] = spk_autocalibration(F,pax);
        end
        if isempty(events)
            % If events is empty, calibration was not possible
            result = struct('quality',nan,'tau',nan,'a',nan,'sigma',nan,'neg',0,'nan',1,'autocal',true);
            return;
        end
        % Replace fixed with calibrated parms
        parms.mlparms.finetune.sigma = sigmaEst;
        parms.mlparms.tau = tauEst;
        parms.mlparms.a = ampEst;
    end


    % Run with parms (or rerun with calibrated parms)
    [spikeTimes,estimatedF,~,parEst] = spk_est(F,parms.mlparms);

    % fit - fit of the F signal at each sample - used to estimate quality
    result.quality = double(corr(estimatedF(~isNaN),F(~isNaN),Type="Pearson"));
    result.tau =parms.mlparms.tau;
    result.a =parms.mlparms.a;
    result.sigma =parEst.finetune.sigma;
    result.neg = mean(isNegative);
    result.nan = mean(isNaN);
    result.autocal = doAutoCal;

    % Save spike times and parameters to file for later use.
    save(filename,"spikeTimes","result");
    assert(exist(filename,"file"),"Could not save mlspike results to  %s",filename);
end
end



