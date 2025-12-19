function [perExpt,perChannel] = artifactDetection(C,parms)
% Several algorithms to detect artifacts in the neural data.
% The output is a list of trials with artifacts. This function is meant to
% be used as the 'fun' in ArtifactParms to fill the Artifact table.
% It is called from populate(ns.Artifact) 
%
% 'Method'
% syncLfp  - detect trials in which the correlation between LFPs on
% is unusually high across a large fraction of electrodes
%   'z': z-scored correlations between electrodes above this
%                       are defined as outliers.
%   'frac' : trials with this fraction of highly correlated
%                       electrodes are outlier trials
%    'time' time window to use
%     'align' align event.
%   'e'  electrodes to use ([] = all loaded)
%
% BK  - Apr 2020
% Updated for ns.C/ns.Artifact
arguments
    C (1,1) ns.C {mustHaveRows(C,1)}
    parms (1,1) struct
end 
% parms 
% .start - signal window start (relative to firstframe)
% .stop  - signal window stop (relative to firstframe, inf is end of trial)

perExptMehods = "sync";
perChannelMethods = ["hasnan" "zscore" "amplitude"];

%% Get the signals per trial
time = fetch1(C,'time');
if numel(time)==3
    step = (time(2)-time(1))/time(3);
else
    step =median(diff(time));
end
%Pull signal at the sampling rate
signal = align(C,removeArtifacts =false,crossTrial =false,start=parms.start,stop=parms.stop,step=step,interpolation="nearest");
signal = timetableToDouble(signal);
trialStartTime = get(ns.Experiment & C, 'cic','prm','firstframe','atTrialTime',inf,'what','clocktime');
trialStopTime= get(ns.Experiment & C,'cic','prm','trialstoptime','what','clocktime');
trialDuration = diff([trialStartTime;trialStopTime(end)])';
[nrSamples,nrTrials,nrChannels] = size(signal);

%% Detect artifacts that apply to all channels
isArtifactTrial = false(1,nrTrials);
for m = intersect(string(fieldnames(parms)),perExptMehods)'
    switch upper(m)
        case 'SYNC'           
            TRIALSNEEDED  =10;
            if nrTrials<TRIALSNEEDED  
                fprintf(2,'Artifact detection based on SYNC needs at least %d trials to get an estimate of normal synchrony.Skipped. \n',TRIALSNEEDED )
                continue;
            end
            nrChannelsInSync = numel(parms.sync.channel);
             r = nan(nrChannelsInSync,nrChannelsInSync,nrTrials);
             for tr=1:nrTrials               
                r(:,:,tr) = corrcoef(squeeze(signal(:,tr,parms.sync.channel)),'rows','pairwise'); % use only trials where each electrode is non-nan
             end             
             r = r -repmat(diag(NaN*ones(1,nrChannelsInSync)),[1 1 nrTrials]);% Remove self corr
             z  =(r-repmat(mean(r,3,"omitnan"),[1 1 nrTrials]))./repmat(std(r,0,3,"omitnan"),[1 1 nrTrials]);
             outlier = z>parms.sync.z;
            nrOutliers = squeeze(sum(outlier,[1 2]))/2; % /2 because correlation is symmetric
            thisIsArtifactTrial = (nrOutliers'>parms.sync.frac*nrChannelsInSync);
        otherwise
            error('Unknown artifact detection method %s\m',m)            
    end
    fprintf('%d experiment artifact trials based on %s\n',sum(thisIsArtifactTrial),m);    
    isArtifactTrial = isArtifactTrial | thisIsArtifactTrial;
end

% These algorithms only identify artifact trials, not time periods.
% Return empty for start/
perExpt.start=[];
perExpt.stop = []; 
perExpt.trial = find(isArtifactTrial);

signal(:,isArtifactTrial,:) =NaN;  % Per-electrode artifact detection will ignore these trials that are already labeled

%% Detect artifacts per channel
% Modes are applied in the order of specification. 
isArtifactTrial = false(nrTrials,nrChannels);
for m = intersect(string(fieldnames(parms)),perChannelMethods)'
    switch upper(m)
        case 'HASNAN'
             out = mean(isnan(signal),1) > parms.hasnan.frac;                  
        case 'ZSCORE'            
             z = (signal-repmat(mean(signal,2,"omitnan"),[1 nrTrials 1]))./repmat(std(signal,0,2,"omitnan"),[1 nrTrials 1]);
             outlier = abs(z)>parms.zscore.z;
             fraction = mean(outlier,1,"omitnan"); % Average over samples 
             out = fraction>parms.zscore.frac;            
        case 'AMPLITUDE'
             isTooHigh= abs(signal)>parms.amplitude.absmax;
             fraction = mean(isTooHigh,1,"omitnan");
             out = fraction>parms.amplitude.frac;            
        otherwise
            error('Unknown artifact detection method %s\m',m)  
    end

    ix = find(reshape(out,nrTrials,nrChannels));
    if isfield(parms.(m),'remove') && parms.(m).remove
        % Remove these artifacts from signal so that they are ignored in
        % the subsequent steps. Defaults to no, user can set with 
        % parms.amplitude.remove = true
        signal(:,ix) =NaN;
    end
    fprintf('%d channel artifact trials found by %s\n',numel(ix),m);    
    isArtifactTrial(ix)= true;
end

[trialNr,channelNr] = find(isArtifactTrial);
 fprintf('%d Artifact trials across %d electrodes\n',numel(unique(trialNr)),numel(unique(channelNr)));    

% These algorithms only identify artifact trials, not time periods.
% Return empty for start/stop
trials  =cellfun(@(x)(find(x)'),num2cell(isArtifactTrial,1),'uni',false);
channelsWithArtifacts = any(isArtifactTrial,1);
perChannel =struct('start',[],'stop',[],'trial',trials(channelsWithArtifacts)','channel',num2cell(find(channelsWithArtifacts),1)');



end