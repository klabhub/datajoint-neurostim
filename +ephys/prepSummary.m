function T = prepSummary(tbl)
% Summarize the original number of channels, the channels that were
% interpolated and the ones that stayed noisy after preprocessing (and were
% discarded by the egi.ephys.read function).
% 
% This applies only to preprocessing using the PREP pipelin in EEGLAB,
% applied to EGI data in ephys.egi.read.
arguments
    tbl (1,1) ns.C    
end

T = fetchtable(tbl,'info');

 isPrep = rowfun(@(x)(isfield(x,'etc')),T,InputVariables="info",OutputFormat="uniform");
if ~all(isPrep)
    fprintf(2,"Some of these ns.C entries do not have an info.etc field:\n");
    T(~isPrep,:);
end

I  = rowfun(@analyzeInfo,T(isPrep,:),InputVariables="info",OutputVariableNames=["nrInterpolated" "nrStillNoisy"],ErrorHandler=@errorFunc);
T= [T(isPrep,:) I];

end


function [nrInterpolated, nrStillNoisy]= analyzeInfo(i)
    nrInterpolated = numel(i.etc.noiseDetection.interpolatedChannelNumbers);
    nrStillNoisy   = numel(i.etc.noiseDetection.stillNoisyChannelNumbers);
end


function [A,B] = errorFunc(me,varargin)
    warning(me.identifier,me.message);
    A = NaN;
    B = NaN;
end