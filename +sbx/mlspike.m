function [varargout] = mlspike(key,cparms,pv)
arguments
    key (1,1) struct
    cparms struct  = struct([]) %#ok<INUSA>
    pv.calibration (1,1) logical = false
end

%assert(exist('bricks.fn_structmerge','file'),"The brick repository must be on the path for mlSpike");
assert(exist('spk_est.m','file'),"The spikes repository must be on the path for mlSpike");
warning('off','backtrace');

parms = fetch(sbx.SpikesParm& struct('stag',key.ctag),'deconv','calibration');
if isfield(parms.deconv,'fluorescence')
    % Identifies a specific ns.C row
    fluorescence = parms.deconv.fluorescence;
else
    % Default name for F
    fluorescence = "fluorescence";
end
% Query the CChannel table for fluorescence time series
allF = (ns.C & struct('ctag',fluorescence) )* (ns.CChannel &  proj(sbx.PreprocessedRoi & key,'roi->channel'));
nrRoi = count(allF);

info = sbx.readInfoFile(fetch(ns.Experiment & key,'LIMIT 1')); % Read one info to get the number of planes
parms.deconv.dt = 1./(fetch1(sbx.Preprocessed & key,'framerate')/info.nrPlanes); % Match dt to framerate
parms.deconv.algo.dogpu = true;
parms.calibration.dogpu =true;
if pv.calibration
    % Call from Mlspikecalibration - run calibration
    assert(all(isfield(parms.calibration,["amin" "amax" "taumin" "taumax" "maxamp" "nrRoi"])),'The %s SpikesParm does not have the required calibration parameters\n',parms.stag)
    nrCalibrationRoi = min(nrRoi,parms.calibration.nrRoi);
    parms.calibration.dt = parms.deconv.dt; % Match to the framerate
    parms.calibration = rmfield(parms.calibration,"nrRoi");
    out = repmat(struct('quality',[],'tau',[],'a',[],'sigma',[]),[nrCalibrationRoi 1]);
    % Get a random subset
    fTpls = fetch(allF,sprintf('ORDER BY rand() LIMIT %d',nrCalibrationRoi));
    pool = nsParPool;
    if ~isempty(pool)
        parfor i=1:nrCalibrationRoi
            out(i) = calibrate(fetch1(ns.CChannel & fTpls(i),'signal'),parms.deconv,parms.calibration);
        end
    else
        for i=1:nrCalibrationRoi           
            out(i) = calibrate(fetch1(ns.CChannel & fTpls(i),'signal'),parms.deconv,parms.calibration);
        end
    end
    varargout{1} = out;
else
    fTpls = fetch(allF);
    nrSamples= fetch1(allF,'nrsamples','LIMIT 1');    
    signal = nan(nrSamples,nrRoi);
    quality = nan(nrRoi,1);
    sigma = nan(nrRoi,1);
    
    if  ~isempty(parms.calibration)
        % This SpikesParm specifies a calibration. Look up the results
        calibration = sbx.Mlspikecalibration & struct('stag',parms.stag) & key;
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
    end
    pool = nsParPool;
    nrRoi =10;
    if ~isempty(pool)
        parfor i=1:nrRoi
            [signal(:,i),quality(i),sigma(i)] = deconvolve(fetch1(ns.CChannel & fTpls(i),'signal'),parms.deconv);
        end
    else
        for i=1:nrRoi            
            F =fetch1(ns.CChannel & fTpls(i),'signal');            
            %F = gpuArray(F); % Tried this but it does not speed up much
           [signal(:,i),quality(i),sigma(i)]  = deconvolve(F,parms.deconv);           
        end
    end
    % Make output for ns.C
    cKey = rmfield(key,"ctag");
    time  =fetch1(ns.C & struct('ctag',fluorescence) & cKey ,'time');  % Same time as the fluorescence
    channelInfo = struct('quality',quality,'sigma',sigma);
    recordingInfo = struct('stag',parms.stag);
    varargout = {signal,time,channelInfo,recordingInfo};
end
end

function [out] = calibrate(F,mlparms, autocalparms)
arguments
    F (:,1) double
    mlparms (1,1) struct
    autocalparms (1,1) struct
end

isNaN = isnan(F);
F(isNaN) = 0; % Could remove samples instead or linearly interpolate, but this should be rare (missing F)

% Autocalibration with the specified parameters
pax = spk_autocalibration('par'); % Get defaults, then overrule with parms.autocalibration
% Copy values from parms.autocalibration struct to pax
pax = fn_structmerge(pax,autocalparms,'strict','recursive','type');
pax.mlspikepar = mlparms;
[tauEst,ampEst,sigmaEst,events] = spk_autocalibration(F,pax);
if isempty(events)
    % If events is empty, calibration was not possible
    out = struct('quality',nan,'tau',nan,'a',nan,'sigma',nan);
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

isNaN = isnan(F);
F(isNaN) = 0; 
[spikeTimes,estimatedF,~,parEst] = spk_est(F,mlparms);
% estimated of the F signal at each sample - used to estimate quality
quality = double(corr(estimatedF(~isNaN),F(~isNaN),Type="Pearson"));
sigma =parEst.finetune.sigma;
signal = sparse(round(spikeTimes./mlparms.dt),1,1,size(F,1),1);
end
