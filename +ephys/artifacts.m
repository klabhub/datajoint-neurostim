function [trials] = artifacts(signal,pk)
% Several algorithms to detect artifacts in the neural data.
%The output is a list of trials with artifacts.
%
% Each method defines an artifact trial differently, and requires different
% structs to specify the definition.  If you specify multilpe structs, each
% will be applied.
%
%  sync : a trial in which the average  absolute value of the correlation 
%               between signals on different electrodes is higher than sync.corr 
%           OR
%           the fraction of electrodes witha a z-score of correlations
%           higher than sync.z is higher than sync.frac
%   
% hasnan: a trial in which more than hasnan.frac fraction of electrodes have a NaN in the specified window. 
%
% clip: An artifact trial in which the abs(signal) of clip.frac electrodes goes above clip.signal 
%           OR
%       the fraction of electrodes witha a z-score of signals
%          higher than clip.z is higher than clip.frac
%   
% BK  - Apr 2020
arguments
    signal timetable
    pk.sync struct  = struct([])
    pk.hasnan  struct = struct([])  
    pk.clip  struct  = struct([])
    pk.start (1,1) = double(min(signal.Time))
    pk.stop (1,1) = double(max(signal.Time))
end

DEFINEDMETHODS = {'sync','hasnan','clip'};

inWindow = timerange(seconds(pk.start),seconds(pk.stop));
signal  = signal(inWindow,:);
[signal,time] = timetableToDouble(signal); %#ok<ASGLU> 
[nrSamples,nrTrials,nrChannels] = size(signal); %#ok<ASGLU> 

methods = find(~cellfun(@(x) isempty(pk.(x)),DEFINEDMETHODS));
nrMethods = numel(methods);
isArtifactTrial = false(nrTrials,1);
for i=1:nrMethods    
    switch upper(DEFINEDMETHODS{methods(i)})
        case 'HASNAN'
             nanInWindow = any(isnan(signal),1);
             fracChannelWithNan = squeeze(mean(nanInWindow,3));
             isArtifactTrial = isArtifactTrial | fracChannelWithNan(:)> pk.hasnan.frac;            
        case 'SYNC'
            signalCorr = nan(nrChannels,nrChannels,nrTrials);
            for tr=1:nrTrials
                signalCorr(:,:,tr) = corrcoef(squeeze(signal(:,tr,:)));
            end
            diagIx = repmat(eye(nrChannels,nrChannels),[1 1 nrTrials]);
            diagIx = logical(diagIx(:));
            signalCorr(diagIx) =NaN; % Remove self corr
            if isfield(pk.sync,'corr')
                meanCorr=squeeze(mean(signalCorr,[1 2],"omitnan"));
                isArtifactTrial = isArtifactTrial |  abs(meanCorr)> pk.sync.corr;
            elseif isfield(pk.sync,'z')
                zLfpCorr  =(signalCorr-repmat(mean(signalCorr,3,"omitnan"),[1 1 nrTrials]))./repmat(std(signalCorr,0,3,"omitnane"),[1 1 nrTrials]);
                outlier = abs(zLfpCorr)>pk.sync.z;
                nrOutliers = squeeze(sum(outlier,[1 2]))/2; % /2 because correlation is symmetric               
                isArtifactTrial = isArtifactTrial |  nrOutliers>pk.sync.frac*nrElectrodes;
            else 
                error('Sync specification error')
            end          
        case 'CLIP'
            if isfield(pk.clip,'signal')
                isArtifactTrial = isArtifactTrial | squeeze(mean(any(abs(signal) > pk.clip.signal,1),3,"omitnan"))'>pk.clip.frac;
            elseif isfield(pk.clip,'z')
                zSignal = (signal-repmat(mean(signal,2,"omitnan"),[1 nrTrials 1]))./repmat(std(signal,0,2,"omitnan"),[1 nrTrials 1]);
                outliers = any(abs(zSignal)>pk.clip.z,1);
                isArtifactTrial = isArtifactTrial |  squeeze(mean(outliers,3))' >pk.clip.frac;         
            else 
                error('Clip specification error')
            end
        otherwise
            error('Unknown methods %s ', DEFINEDMETHODS{i})
    end    
end
trials =find(isArtifactTrial);
end