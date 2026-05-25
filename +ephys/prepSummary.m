function T = prepSummary(tbl,pv)
% Summarize the original number of channels, the channels that were
% interpolated and the ones that stayed noisy after preprocessing (and were
% discarded by the egi.ephys.read function).
% 
% This applies only to preprocessing using the PREP pipelin in EEGLAB,
% applied to EGI data in ephys.egi.read.
arguments
    tbl (1,1) ns.C    
    pv.plot (1,1) logical = false
    pv.nrTotalExpected (1,1) double = 257
end

T = fetchtable(tbl,'info');
Nch = fetchtable(aggr(tbl, ns.CChannel, 'count(channel)->nrChannels'), 'nrChannels');

 isPrep = rowfun(@(x)(isfield(x,'etc')),T,InputVariables="info",OutputFormat="uniform");
if ~all(isPrep)
    fprintf(2,"Some of these ns.C entries do not have an info.etc field:\n");
    T(~isPrep,:);
end

I  = rowfun(@analyzeInfo,T(isPrep,:),InputVariables="info",OutputVariableNames=["nrInterpolated" "nrStillNoisy"],ErrorHandler=@errorFunc);
T= [T(isPrep,:) I];
T = innerjoin(T,Nch);
T =convertvars(T,@isnumeric,"double"); % Get rid of int64

nrTotal = T.nrStillNoisy + T.nrChannels;
unexpected = nrTotal ~=pv.nrTotalExpected;
if any(unexpected)
    fprintf(2, 'Unexpected channel numbers:\n')
    disp(T(unexpected,:))
end

if pv.plot
    
    figByName('PREP Summary')
    clf;
    nexttile
    percentage = 100*T.nrInterpolated(~unexpected)./nrTotal(~unexpected);
    histogram(percentage,0:5:100)
    xlabel '%Interpolated'
    ylabel '#Experiments'
    title (sprintf('Median: %.1f%% IQR: %.1f%%',median(percentage),iqr(percentage)))
    
    nexttile
    percentage =100*T.nrStillNoisy(~unexpected)./nrTotal(~unexpected);
    histogram(percentage)
    xlabel '%Removed (still noisy)'
    ylabel '#Experiments'
    title (sprintf('Median: %.1f%% IQR:%.1f%%',median(percentage),iqr(percentage)))   
end

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