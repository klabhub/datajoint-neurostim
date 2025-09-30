function [varargout] = mlspike(key,cparms,pv)
arguments
    key (1,1) struct
    cparms (1,:) struct  = struct([]) %#ok<INUSA> % This is not used, just a placeholder so that it can be called from ns.C/maketuples
    pv.calibration (1,1) logical = false
end

assert(exist('xplor.m','file'),"The xplor repository must be on the path for sbx.mlspike");
assert(exist('spk_est.m','file'),"The spikes repository must be on the path for sbx.mlspike");
warning('off','backtrace');

parms = fetch(sbx.SpikesParm& struct('stag',key.ctag),'deconv','calibration','fluorescence');
% Query the CChannel table for fluorescence time series
allF = (ns.C & struct('ctag',parms.fluorescence) )* (ns.CChannel &  proj(sbx.PreprocessedRoi & key,'roi->channel'));
nrRoi = count(allF);
prep = fetch(sbx.Preprocessed & key,'framerate','nrplanes');
parms.deconv.dt = 1./(prep.framerate/prep.nrplanes); % Match dt to framerate
nrSamples= fetch1(allF,'nrsamples','LIMIT 1');

if pv.calibration
    % Call from Mlspikecalibration - run calibration
    assert(all(isfield(parms.calibration,["amin" "amax" "taumin" "taumax" "maxamp" "nrRoi"])),'The %s SpikesParm does not have the required calibration parameters\n',parms.stag)
    nrRoi = min(nrRoi,parms.calibration.nrRoi);
    parms.calibration.dt = parms.deconv.dt; % Match to the framerate
    parms.calibration = rmfield(parms.calibration,"nrRoi");
    % Get a random subset of ROIS and prepare arguments for the parfeval
    fTpls = fetch(allF,sprintf('ORDER BY rand() LIMIT %d',nrRoi));
    fun= @calibrate;
    nrOut = 1;
    inputArgs = {parms.deconv,parms.calibration};
else
    fTpls = fetch(allF);   
    if  ~isempty(parms.calibration)
        % This SpikesParm specifies a calibration. Look up the results
        calibration = sbx.Mlspikecalibration & struct('stag',parms.stag) & key;
       %{
        assert(exists(calibration),"No %s calibration found for %s on %s at %s. \n Populate the sbx.Mlspikecalibration table first. ",parms.stag,key.subject,key.session_date,key.starttime);
        calibration =fetch(calibration,'*');
        if isnan(calibration.quality)
            % All calibration failed, but at least we tried.  Use the
            % default parameters as speci
        else
            parms.deconv.a = calibration.a;
            parms.deconv.tau = calibration.tau;
            parms.deconv.finetune.sigma = calibration.sigma;
        end
       %}
    end
    fun= @deconvolve;
    nrOut = 3;
    inputArgs = {parms.deconv};
end

%Do the work on a parpool (or serial if nsParPool returns empty)
pool = nsParPool;
fprintf('Queue %d channels with %d samples for %s at %s\n',nrRoi,nrSamples,func2str(fun),datetime("now"))
for i=1:5 %nrRoi
    future(i) = parfeval(pool,fun,nrOut,fetch1(ns.CChannel & fTpls(i),'signal'),inputArgs{:}); %#ok<AGROW>
end
afterEach(future,@done,0,PassFuture = true);  % Update the command line
wait(future); % Wait for all to complete
keep = cellfun(@isempty,{future.Error});% Remove errors
assert(any(keep),'All %s failed on %s on %s @ %s ',func2str(fun),key.subject,key.session_date);
if pv.calibration
    % Single struct output
    varargout{1} = fetchOutputs(future(keep));    % Collect results
else
    % Deconvolve call from ns.C/maketuples 
    [signal,quality,sigma] = fetchOutputs(future(keep),UniformOutput=false);    % Collect results
    % Convert to necessary output for ns.C/makeTuples
    cKey = rmfield(key,"ctag");
    time  =fetch1(ns.C & struct('ctag',parms.fluorescence) & cKey ,'time');  % Same time as the fluorescence
    channelInfo = struct('quality',quality,'sigma',sigma);
    recordingInfo = struct('stag',parms.stag);
    varargout = {cat(2,signal{:}),time,channelInfo,recordingInfo};
end


warning('on','backtrace');

    function done(ftr)
        % Message to comman line that a job is done or errored out
        if ~isempty(ftr.Error)
            fprintf("%s failed after %s - (Error: %s)\n ",func2str(ftr.Function),ftr.RunningDuration,ftr.Error);
            ftr.Diary
        else
            fprintf("%s is done after %s - (State: %s)\n ",func2str(ftr.Function),ftr.RunningDuration,ftr.State);
        end
    end
end

function [out] = calibrate(F,mlparms, autocalparms)
arguments
    F (:,1) double
    mlparms (1,1) struct
    autocalparms (1,1) struct
end
global BRICKPROGRESS %#ok<GVMIS>
BRICKPROGRESS = false;

isNaN = isnan(F);
F(isNaN) = 0; % Could remove samples instead or linearly interpolate, but this should be rare (missing F)
isNegative  = F<0;
% F is supposed to be F/F0 and therefore positive, but the neuropil
% correction can occasionally generate negative F.  Set those to zero.
if any(isNegative)
    fprintf('%.2f %% of samples have negative F. Set to 0.\n',100*mean(isNegative));
    F(isNegative) = 0;
end

% Autocalibration with the specified parameters
pax = spk_autocalibration('par'); % Get defaults, then overrule with parms.autocalibration
% Copy values from parms.autocalibration struct to pax
pax = brick.structmerge(pax,autocalparms,'strict','recursive','type');
pax.mlspikepar = mlparms;
pax.mlspikepar.dographsummary = false;
[tauEst,ampEst,sigmaEst,events] = spk_autocalibration(F,pax);
if isempty(events)
    % If events is empty, calibration was not possible
    out = struct('quality',nan,'tau',nan,'a',nan,'sigma',nan,'neg',0,'nan',1);
else
    % Rerun with calibrated parms
    pax = mlparms;
    % Replace fixed with calibrated parms
    pax.finetune.sigma = sigmaEst; % Always estimated
    pax.tau = tauEst;
    pax.a = ampEst;
    [~,estimatedF,~,parEst] = spk_est(F,pax);
    % fit - fit of the F signal at each sample - used to estimate quality
    out.quality = double(corr(estimatedF(~isNaN),F(~isNaN),Type="Pearson"));
    out.tau =tauEst;
    out.a =ampEst;
    out.sigma =parEst.finetune.sigma;
    out.neg = mean(isNegative);
    out.nan = mean(isNaN);
end
end

function [signal,quality,sigma] = deconvolve(F,mlparms)
% This function performs spike deconvolution using the Deneux et al mlspike algorithm,
% tpl  -  ns.CChannel table row that containst the (neuropil corrected)  fluorescence data.
% mlparms - deconvolution model parameters
% autocalparms - autocalibration parameters
%
% Note that MLspikecalibration/makeTuples calls this with autocalparms and
% mlparms from Mlspikeparm table (to do the calibration), while  Spikes/makeTuples
% calls it with just mlparms (modified by the calibrated values)
%
%
arguments
    F (:,1) double
    mlparms (1,1) struct   
end
global BRICKPROGRESS %#ok<GVMIS>
BRICKPROGRESS = false;

isNaN = isnan(F);
F(isNaN) = 0;
isNegative  = F<0;
% F is supposed to be F/F0 and therefore positive, but the neuropil
% correction can occasionally generate negative F.  Set those to zero.
if any(isNegative)
    fprintf('%.2f %% of samples have negative F. Set to 0.\n',100*mean(isNegative));
    F(isNegative) = 0;
end
[spikeTimes,estimatedF,~,parEst] = spk_est(F,mlparms);
% estimated of the F signal at each sample - used to estimate quality
quality = double(corr(estimatedF(~isNaN),F(~isNaN),Type="Pearson"));
sigma =parEst.finetune.sigma;
signal = sparse(round(spikeTimes./mlparms.dt),1,1,size(F,1),1);
end
