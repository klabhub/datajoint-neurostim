function [signal,time,channelInfo,recordingInfo] = mlspike(key,parms)
% Function that uses the spikes and xplor repositories from Deneux et al.
% 2016 to run spike deconvolution.
% 
% This function is called when populating ns.C with the appropriate
% parameters. 
%
% With autocalibration defined, one calibration is run (with nrRoiForCal
% randomly chosen rois), and the tau, amplitude, and sigma that result from
% this are used for all subsequent roi. This is done to reduce the time spent
% on calibration. If you prefer to run autocalibration on every roi, set
% nrRoiForCal to inf. 
% 
% EXAMPLE
% Setup the parameters for the alogrithm (see spikes/spk_demo.m)
%
%mlparms = sbx.mlspikeDefaults("DENEUX16","gcamp6s");
%mlparms.algo.nspikemax =4;  % Allow 4 spikes per bin (15 Hz)
%mlparms.dographsummary = false;
%mlparms.display = 'none';
% 
%autocalparms.amin = 0.05;
%autocalparms.amax =0.2;
%autocalparms.taumin = 0.25;                 
%autocalparms.taumax = 2;
%autocalparms.maxamp =4  % A maxamp of 4 seems necessary in our data.
%autocalparms.display = 'none'; % No figures
%
% Create a ns.CParm struct that uses these parms
% deneux = struct('ctag','deneux16',...       % Name of this C 
%                   'description','Deneux spikes',... 
%                    'extension','.sbx',...  
%                    'fun','sbx.mlspike',...      
%                    'parms',struct( 'prep','gcamp6s',... % Which preprocessed set to use (links to sbx.Preprocessed)
%                                    'restrict','pcell>0.75 AND radius>2.185',...
%                                    'mlparms',mlparms,...
%                                    'autocal',autocalparms,...
%                                     'nrRoiForCal',20)); % Use 20 ROIs per session to calibrate 
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
    filename = fullfile(fldr,"plane" + string(pl), string(roi) + "." + key.ctag + ".mlspike.mat");
    done = false(nrRoiThisPlane,1);
    for i=1:nrRoiThisPlane
        done(i) = exist(filename(i),"file");
    end
    plane = repmat(pl,nrRoiThisPlane,1);
    roiToDo = [roiToDo ; table(filename,done,roi,plane)]; %#ok<AGROW>
end
if any(~roiToDo.done)
    % For the files that do not yet exist, check whether a calibration exists
    % (if requested)
    needCalibration = isfield(parms,'autocal') && ~isempty(parms.autocal);
    if needCalibration
        haveCalibration = exists(sbx.Mlcalibration & key);
        if ~haveCalibration
            % Run calibration for the session
            nrRoi = height(roiToDo);
            assert(all(isfield(parms.autocal,["amin" "amax" "taumin" "taumax" "maxamp"])),'The %s does not have the required calibration parameters\n',key.ctag)
            if isfield(parms,'nrRoiForCal')
                nrRoi = min(nrRoi,parms.nrRoiForCal);
            end
            % Select rois from the roiT
            useRoi = randperm(height(roiToDo),nrRoi);
            % Calibrate based on this subset
            calResults = loop(roiToDo(useRoi,:),parms,fldr,true);
            % This only has the parameters (and the quality); the data have
            % been saved to disk. 
            % Store the calibration result
            tpl = ns.stripToPrimary(sbx.Preprocessed,prep);
            tpl.ctag= key.ctag;
            tpl.tau = mean([calResults.tau],"omitmissing");
            tpl.sigma  =mean([calResults.sigma],"omitmissing");
            tpl.a = mean([calResults.a],"omitmissing");
            tpl.quality = mean([calResults.quality],"omitmissing");
            tpl.failed = sum(isnan([calResults.quality]));
            tpl.neg = mean([calResults.neg],"omitmissing");
            tpl.nan = mean([calResults.nan],"omitmissing");
            insert(sbx.Mlcalibration, tpl);
        end
        % Calibration should be available now
        calibration = fetch(sbx.Mlspikecalibration & key, '*');
        assert(~isempty(calibration),"ML Spike calibration must have failed?");
        % Replace values in mlparms with calibrated values
        parms.mlparms.tau = calibration.tau;
        parms.mlparms.a  =calibration.a;
        % parms.mlparms.finetune.sigma = Not set- reestimated for each ROI (cheap);
    else
        % Run uncalibrated. Keep parms.mlparms as is.
    end
    % Recheck files as the calibration may have created some
    for i=1:height(roiToDo)
        roiToDo.done(i) = exist(roiToDo.filename(i),"file");
    end
    deconResults = loop(roiToDo(~roiToDo.done,:),parms,fldr,false); 
    mQuality = mean([deconResults.quality],"omitmissing");
    fprintf('Completed mlspike deconvolution on %d rois. Mean quality %.2f\n',height(roiToDo),mQuality)
end

%% Check that all are now done
for i=1:height(roiToDo)
    roiToDo.done(i) = exist(roiToDo.filename(i),"file");
end
if any(~roiToDo.done)
    fprintf('The following ROIs were not deconvolved:\n')
    roiToDo(~roiToDo.done,:)
end
roiToDo = roiToDo(roiToDo.done,:); % Crop
%% Experiment specific
% Now files with spiketimes exist on disk. Load and extract the relevant frames for the
% current experiment to return to the ns.C caller.

% Use the metadata added to the experiment by sbx.Preprocessed to match
% frames to experiments
exptT = ns.getMeta(ns.Experiment & key ,["nrframes" "nrplanes"]);
exptT = sortrows(exptT,"starttime");
exptT = convertvars(exptT,["nrframes" "nrplanes"],"double");
row = find(key.starttime==exptT.starttime);
cumFrames = cumsum(exptT.nrframes);
if row>1
    start = cumFrames(row-1);
else
    start = 0;
end
keepFrameIx =start + (1:exptT.nrframes(row));
frameNsTime = get(ns.Experiment & key,'mdaq','prm','laserOnDigHigh','what',"clocktime");
nrTTL = numel(frameNsTime);
% There always appears to be 1 extraneous TTL; check the match with
% this assumption
assert((exptT.nrframes(row)+1)==floor(nrTTL/exptT.nrplanes(row)),'Cannot map SBX frames to trials; TTL-Frame mismatch (%d TTL %d frames in sbx).\n',nrTTL,exptT.nrframes(row));        
frameNsTime(1) =[];      
signal = [];
nrFrames = numel(frameNsTime);

channelInfo =  struct('nr',num2cell(roiToDo.roi),'quality',NaN,'tau',NaN,'sigma',NaN,'neg',NaN,'nan',NaN);
% Note that ROI are numbered across planes. This is matched in
% sbx.PreprocessedRoi to allow inner joins with CChannel.
% (see example in sxb.PreprocessedRoi/plotSpatial)
for i = 1:height(roiToDo)
    %% Read saved results
    fprintf('Reading spikeml mat files...\n')      
    load(roiToDo.filename(i),'spikeTimes','result');
    channelInfo(i).quality = result.quality;
    channelInfo(i).tau= result.tau;
    channelInfo(i).a= result.a;
    channelInfo(i).neg= result.neg;
    channelInfo(i).nan = result.nan;
    % Convert spike times to a signal with counts
    thisSignal = brick.timevector(spikeTimes,(0:nrFrames-1)*parms.mlparms.dt,'count');
    framesNotInFile = sum(keepFrameIx > numel(thisSignal));
    assert(framesNotInFile==0,"%s has %d too few frames for %s on %s",roiToDo.filename(i),framesNotInFile,key.starttime, key.session_date);
    thisSignal = thisSignal(keepFrameIx);        
    signal = [signal  thisSignal]; %#ok<AGROW>

    fprintf('Done in %s.\n',seconds(toc))
end

time = [frameNsTime(1) frameNsTime(end) nrFrames];
recordingInfo = struct('dummy',true);

end

function [out] = loop(roiT,parms,fldr,calibrate)
arguments
    roiT table % List of rois that have not been done yet
    parms (1,1) struct
    fldr (1,1) string
    calibrate (1,1) logical
end
global BRICKPROGRESS %#ok<GVMIS>
BRICKPROGRESS = false;

F = [];
for pl = unique(roiT.plane)
    thisFile = fullfile(fldr,"plane" +  string(pl) ,'F.npy');
    if ~exist(thisFile,"file")
        error('File %s does not exist',thisFile);
    end
     thisSignal =  ndarrayToArray(py.numpy.load(thisFile,allow_pickle=true),single=true);
    keepRoi = roiT.roi(roiT.plane==pl);
    F = [F;thisSignal(keepRoi,:)']; %#ok<AGROW>
end
%% Use parallel pool if requested
pool = nsParPool;
[nrSamples, nrRoi] = size(F);
fprintf('Queue %d channels with %d samples for deconvolution at %s\n',nrRoi,nrSamples,datetime("now"))
tStart = tic;
nrRoiDone = 0;

for i=1:nrRoi
    future(i) = parfeval(pool,@deconvolve,1,F(:,i),parms,roiT.filename(i),calibrate); %#ok<AGROW>
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

function result = deconvolve(F,parms,filename,calibrate)
% Function that does the actual mlspike deconvolution.
arguments
    F (:,1) single  % Fluorescence signal
    parms (1,1) struct   % parms
    filename  (1,1) string
    calibrate (1,1) logical = false  % Set to true to perform autocal 
end


isNaN = isnan(F);
F(isNaN) = 0; % Could remove samples instead or linearly interpolate, but this should be rare (missing F)
isNegative  = F<0;
% F is supposed to be F/F0 and therefore positive, but the neuropil
% correction can occasionally generate negative F.  Set those to zero.

if mean(isNegative)>0.05
    % Fail if more than 5% are negative
    result = struct('quality',0,'tau',nan,'a',nan,'sigma',nan,'neg',mean(isNegative),'nan',mean(isNaN));
    return;
else
    if any(isNegative)
        fprintf('%.2f %% of samples have negative F. Set to 0.\n',100*mean(isNegative));
        F(isNegative) = 0;
    end

    if calibrate
        % Autocalibration with the specified parameters
        pax = spk_autocalibration('par'); % Get defaults, then overrule with parms.autocal
        % Copy values from parms.autocal struct to pax
        pax = brick.structmerge(pax,parms.autocal,'strict','recursive','type');        
        pax.dt = parms.mlparms.dt;
        pax.mlspikepar = parms.mlparms;
        pax.mlspikepar.dographsummary = false;
        [tauEst,ampEst,sigmaEst,events] = spk_autocalibration(F,pax);
        if isempty(events)
            % If events is empty, calibration was not possible
            result = struct('quality',nan,'tau',nan,'a',nan,'sigma',nan,'neg',0,'nan',1);
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

    % Save spike times and parameters to file for later use.
    save(filename,"spikeTimes","result");
    assert(exist(filename,"file"),"Could not save mlspike results to  %s",filename);

end
end



