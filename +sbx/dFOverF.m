function [dFF,shotNoise]  = dFOverF(tFluorescence,tFNeuropil,baselineWindow,pv)
arguments
    tFluorescence
    tFNeuropil
    baselineWindow
    pv.percentile = 8
end

tCorrected = tFluorescence;
if ~isempty(tFNeuropil)
    tCorrected.Variables  = tCorrected.Variables - 0.7* tFNeuropil.Variables; % Subtract neuropil
end

[nrTimePoints, nrTrials] = size(tCorrected);
nrRoi = numel(tCorrected{1,1});

F = baseline(tCorrected,baselineWindow,pv.percentile);
F = repmat(F,[nrTimePoints 1]);
dF  = table2array(tCorrected);
dF  = reshape(dF,[nrTimePoints nrRoi nrTrials]);
dFF = 100*(mean(dF,3,"omitnan")-F)./F;

%% Shot noise
% Rupprecht et al. A measure of shot noise.
if nargout>1
    % Compute df/f without neuropil subtraction
    fUncorrected = baseline(tFluorescence,baselineWindow,pv.percentile);
    fUncorrected = repmat(fUncorrected,[nrTimePoints 1 nrTrials]);
    dFUncorrected  = table2array(tFluorescence);
    dFUncorrected    = reshape(dFUncorrected  ,[nrTimePoints nrRoi nrTrials]);
    dFFUncorrected   = 100*(dFUncorrected-fUncorrected)./fUncorrected;
    % Then compare subsequent frames and scale by the sqrt of the
    % framerate (averaging over trials and time points : [1 3])
    shotNoise  = median(abs(diff(dFFUncorrected,1)),[1 3],"omitnan")/sqrt(1./seconds(median(diff(tCorrected.Time))));
end

end

function F = baseline(T,baselineWindow,percentile)
%% Compute baseline as a percentile of the fluorescnece in a time window
inBaseline = timerange(seconds(baselineWindow(1)),seconds(baselineWindow(2)));
tBaseline   = T(inBaseline,:);
[nrTimePointsBaseline,nrTrials] =size(tBaseline);
nrRoi = numel(T{1,1});
F = table2array(tBaseline);
F = reshape(F,[nrTimePointsBaseline nrRoi nrTrials]);
F = prctile(F,percentile,[1 3]);  % Determine percentile over time points in the window and all trials
end
