function [dFF,shotNoise]  = dFOverF(tFluorescence,tFNeuropil,baselineWindow,pv)
% Compute dF/F from fluorescence and  neuropil
% INPUT
% tFluorescence - Timetable with fluorescence data. [nrTimePoints nrTrials]
% with each entry representing the fluorescence of each ROI [1 nrRoi].
% tFNeuropil - Timetable with neuropil fluorescence for each ROI. 
%  baselinwWindow - Start and Stop of the window that defines the baseline.
% PArm/Value Pairs:
% percentile  - Percentile of the histogram of fluorescence values that is
% considered the background (F0). [8].
% neuropilFactor - Neuropil correction uses 
%                   F ->F-neuropilFactor*Fneuropil.  [0.7]
% OUTPUT
% dFF - dF/F  [nrTimePoints nrRoi]
% shotNoise - Shot Noise estimate per ROI.
% 
% NOTES
%  Fluorescence is first corrected for Neuropil F (if not empty), then the
%  baseline (F0) is determined as the pv.percentile of the dstribution of
%  corrected F in the baseline window. dF/F is then the mean of the
%  corrected fluorescence acorss trials minus the F0 and then divided by
%  F0.

arguments
    tFluorescence
    tFNeuropil
    baselineWindow
    pv.percentile = 8
    pv.neuropilFactor = 0.7
end

tCorrected = tFluorescence;
if ~isempty(tFNeuropil)
    tCorrected.Variables  = tCorrected.Variables - pv.neuropilFactor* tFNeuropil.Variables; % Subtract neuropil
end

[nrTimePoints, nrTrials] = size(tCorrected);
nrRoi = numel(tCorrected{1,1});

F0 = baseline(tCorrected,baselineWindow,pv.percentile);
F0 = repmat(F0,[nrTimePoints 1 nrTrials]);
dF  = table2array(tCorrected);
dF  = reshape(dF,[nrTimePoints nrRoi nrTrials]);
dFF = mean(100*(dF-F0)./F0,3,"omitnan");


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
    shotNoise  = median(abs(diff(dFFUncorrected,1)),[1 3],"omitnan")/sqrt(1./seconds(median(diff(tFluorescence.Time))));
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
if percentile ==0
    % Used the mean over time points and trials
   F = mean(F,[1 3],"omitnan");   
else
    % Determine mean over trials, and then the percentile over the time
    % points
    %F = prctile(mean(F,3,"omitnan"),percentile);  
    F = prctile(F,percentile,[1 3]);   % Percentile over time and trial

end
end
