function [out] = mlspikecompute(tpl,mlparms,autocalparms)
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
    tpl (1,1) struct 
    mlparms (1,1) struct
    autocalparms (1,1) struct = struct([]);
end

assert(exist('fn_structmerge','file'),"The brick repository must be on the path for mlSpike");
assert(exist('spk_est.m','file'),"The spikes repository must be on the path for mlSpike");
setupPython;
warning('off','backtrace');

    tic;
    F = tpl.signal;
    isNaN = isnan(F);
    F(isNaN) = 0; % Could remove samples instead or linearly interpolate, but this should be rare (missing F)
   
    % perform auto-calibration
    tic
    if ~isempty(autocalparms)
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
    else
        % No calibration - use parms to deconvolve

    end 
end