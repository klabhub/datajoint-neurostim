function [signal,time,channelInfo,recordingInfo] = mlspike(key,parms)
% The key refers to an sbx file in an experiment.
%
%  This function performs spike deconvolution using the Deneux  mlspike algorith,
%
% It first runs autocal on a subset of the ROIs,then uses the median autocal values
% to the rest.

% The parms struct contains the following fields
% .prep -  A unique name to identify the preprocessing instructions (= a row in sbx.PreprocessedParm)
% .fluorescence - e.g. 'fluorescence'    - The ctag that identifies the Fluoresence data
% (neuropil corrected) in the ns.C table that will form the basis for the
% deconvolution.
% .restrict - A restriction on sbx.PreprocessedRoi to limit the ROI to a subset.
%               For instance 'pcell>0.75'

preprocessed = sbx.Preprocessed & key & struct('prep', parms.prep);
assert(exists(preprocessed),'No preprocessed data for %s in session %s for subject %s. Run populate(sbx.Preprocessed,prep="%s") first',parms.prep,key.session_date,key.subject,parms.prep);
fKey = key;
fKey.ctag = parms.fluorescence;
assert(exists(ns.C & fKey ),'No fluoresence data for %s in session %s for subject %s. Run populate(ns.C,ctag="%s") first',parms.prep,key.session_date,key.subject,parms.fluorescence);
assert(exist('fn_structmerge','file'),"The brick repository must be on the path for mlSpike");
assert(exist('spk_est.m','file'),"The spikes repository must be on the path for mlSpike");

setupPython;
warning('off','backtrace');

%% Autocalibration
roi = sbx.PreprocessedRoi & key & parms.restrict;
nrRoi = count(roi);
nrCalibrationRoi = round(max(0.05*nrRoi,min(nrRoi,10)));
calibrationChannels = fetchn(roi,'roi');
calibrationChannels = calibrationChannels(randperm(nrRoi));
calibrationChannels = calibrationChannels(1:nrCalibrationRoi);
fluorescence = (ns.C& struct('ctag',parms.fluorescence) ) *(ns.CChannel & proj(roi,'roi->channel'));

dt = 1/fetch1(preprocessed,'framerate');
% Setup the parameters for autocalibration
pax = spk_autocalibration('par'); % Get defaults, then overrule with parms.autocalibration
% Copy values from parms.autocalibration struct to pax
pax = fn_structmerge(pax,parms.autocalibration,'strict','recursive','type');
% For consistency of autocalibration and deconvolution;
% copy deconv to mlspike par, the autocalibration function
% calls tps_mlspikes with these parms
if isfield(parms,'mlspikepar')
    pax.mlspikepar = parms.mlspikepar;
elseif isfield(parms,'publication')
    pax.mlspikepar = defaults(parms.publication,parms.indicator,parms.finetuned);
else
    error('');
end
pax.dt = dt; 



calCntr = 0;
sigma = nan(nrCalibrationRoi,1);
tau = nan(nrCalibrationRoi,1);
amp = nan(nrCalibrationRoi,1);
quality = nan(nrCalibrationRoi,1);
for  ch = calibrationChannels(:)'
    calCntr =calCntr+1;
    tic;
    F  = fetch1(fluorescence & struct('channel',ch),'signal');
    F = F(1:10000);
    isNaN = isnan(F);
    F(isNaN) = 0; % Could remove samples instead or linearly interpolate, but this should be rare (missing F)
   
    % perform auto-calibration
    tic
    [tauEst,ampEst,sigmaEst,events] = spk_autocalibration(F,pax);
    toc
    % If events is empty, calibration was not possible, skip.
    if ~isempty(events)
        pax.finetune.sigma = sigmaEst; % Always estimated
        pax.tau = tauEst;
        pax.a = ampEst;

        % Do the deconvolution
        tic
        [~,estimatedF,~,parEst] = spk_est(F,pax);
toc
        % fit - fit of the F signal at each sample - used to estimate quality
        quality(calCntr) = double(corr(estimatedF(~isNaN),F(~isNaN),Type="Pearson"));
        tau(calCntr) =tauEst;
        amp(calCntr) =ampEst;
        sigma(calCntr) =parEst.finetune.sigma;

    end
end

%% Find the values to use for all future deconv.


%    nrFrames = size(signal,1);
%    time = [frameNsTime(1) frameNsTime(end) nrFrames];
% Note that ROI are numbered across planes. This is matched in
% sbx.PreprocessedRoi to allow inner joins with CChannel.
% (see example in sxb.PreprocessedRoi/plotSpatial)
channelInfo =  struct('nr',num2cell(rois)');
end


% A function to setup default parameters
function mlspikepar = defaults(publication,indicator,fineTuned)
arguments
    publication (1,1) string
    indicator (1,1) string
    fineTuned (1,1) logical =false
end

mlspikepar = tps_mlspikes('par');
mlspikepar.algo.nspikemax = 4;
switch (publication)
    case "RUPPRECHT25"
        % Return parameters for Gcamp6s based on Rupprecht et al. 2025.
        switch upper(indicator)
            case "GCAMP6S"
                mlspikepar.hill = 1.84;
                mlspikepar.a = 0.113;  % Amplitude
                mlspikepar.tau  =1.87; % Decay tau in s
                mlspikepar.ton = 0.07; % Rise tau (t_on) in s
                mlspikepar.pnonlin = 0.1;
                mlspikepar.drift.parameter = 0.1;
            case "GCAMP8S"
                mlspikepar.hill = 2.2;
                mlspikepar.a = 0.576;  % Amplitude
                mlspikepar.tau  =0.267; % Decay tau in s
                mlspikepar.ton = 0.00472; % Rise tau (t_on) in s
                mlspikepar.pnonlin = 0.1;
                mlspikepar.drift.parameter = 0.1;
                if fineTuned
                    % Use tau and amplitude as estimated by Rupprecht
                    % with a grid search on ground truth data
                    mlspikepar.a    = 1;  % Amplitude
                    mlspikepar.tau  = 0.7; % Decay tau in s
                end
            otherwise
                error("No %s defaults for %s",indicator,publication)
        end
    case "DENEUX16"
        switch upper(indicator)
            case "GCAMP6S"
                mlspikepar.hill = 1; % Using polynomial pnonlin instead
                mlspikepar.a = 0.07;  % Amplitude
                mlspikepar.tau  =1.3; % Decay tau in s
                mlspikepar.ton = 0.02; % Rise tau (t_on) in s
                mlspikepar.pnonlin =  [0.73, -0.05];  % Deneux values for Gcamp6s;
                mlspikepar.drift.parameter = 0.1;
            otherwise
                error("No %s defaults for %s",indicator,publication)
        end
    otherwise
        error("Unknown publication %s",publication)
end
end
