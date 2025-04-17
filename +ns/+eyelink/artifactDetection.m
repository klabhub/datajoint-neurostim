function [perExpt,perChannel] = artifactDetection(C,parms)
% Several algorithms to detect artifacts in eyelink data.
% The output is a list of trials with artifacts. This function is meant to
% be used as the 'fun' in ArtifactParms to fill the Artifact table.
% It is called from populate(ns.Artifact)
%
%
% BK  - Dec  2023

arguments
    C (1,1) ns.C {mustHaveRows(C,1)}
    parms (1,1) struct
end
% parms
% .start - signal window start (relative to firstframe)
% .stop  - signal window stop (relative to firstframe, inf is end of trial)
% .speed.max
% .speed.frac
% .acceleration.max
% .acceleration.frac
% .pupil.frac

perExptMethods = ["SPEED" "ACCELERATION"];
nrTrials =fetch1(ns.Experiment &C,'nrtrials');
isArtifactTrial = false(1,nrTrials);

definedMethods = upper(string(fieldnames(parms))); 
    
if any(ismember(definedMethods,perExptMethods))
    %% Get the x and y position per trial at the sampling rate
    [xy,~] = align(C,channel = {'px','py'}, removeArtifacts =false,crossTrial =false,start=parms.start,stop=parms.stop,interpolation="nearest");
    [nrSamples,nrTrials,nrChannels] = size(xy); %#ok<ASGLU>
    r = sqrt(sum(xy.^2,3));
    for m = definedMethods
        switch upper(m)
            case 'SPEED'
                % TOO FAST
                s = abs(diff(r,1,1));
                outlier = s>parms.speed.max;
                fraction = squeeze(mean(outlier,1,"omitnan")); % Average over samples
                thisIsArtifactTrial = (fraction>parms.speed.frac);
            case 'ACCELERATION'
                % Too much acceleration
                a =diff(r,2,1);
                outlier = a>parms.acceleration.max;
                fraction = squeeze(mean(outlier,1,"omitnan")); % Average over samples
                thisIsArtifactTrial = (fraction>parms.acceleration.frac);
            otherwise
                error('Unknown artifact detection method %s\m',m)
        end
        fprintf('%d experiment artifact trials based on %s\n',sum(thisIsArtifactTrial),m);
        isArtifactTrial = isArtifactTrial | thisIsArtifactTrial;
    end
end

perChannel = struct('channel',[],'start',[],'stop',[]);
%% Pupil based
if ismember("PUPIL",definedMethods) 
    % Uses:
    % pupil.frac = max fraction of samples without change in pupil
    %           area (i.e. no signal)
    % pupil.blink.pre  - samples before startblink event to set to NaN
    % pupil.blink.post - samples after endblink event to set to NaN
    %
    if isfield(parms.pupil,'blink')
        blinkStartTime = get(ns.Experiment & C,'edf','prm','startblink','what','clocktime');        
        perChannel.channel  = fetch1(ns.CChannel & C & 'name="pa"','channel');          
        if ~isempty(blinkStartTime)
            blinkStopTime = get(ns.Experiment & C,'edf','prm','endblink','what','clocktime');
            perChannel.start = [perChannel.start; blinkStartTime - parms.pupil.blink.pre];
            perChannel.stop = [perChannel.stop ;blinkStopTime + parms.pupil.blink.post];
        end       
    end
    if isfield(parms.pupil,'frac')
        % Remove complete trials if a large frac has a pupil area that does
        % not change
        [area,~] = align(C,channel = 'pa', removeArtifacts =false,crossTrial =false,start=parms.start,stop=parms.stop,interpolation="nearest");
        [nrSamples,nrTrials,nrChannels] = size(area); %#ok<ASGLU>
        noChange = diff(area)==0;
        fraction = mean(noChange,1);
        thisIsArtifactTrial = (fraction>parms.pupil.frac);
        fprintf('%d experiment artifact trials based on %s\n',sum(thisIsArtifactTrial),'Pupil');
        isArtifactTrial = isArtifactTrial | thisIsArtifactTrial;
    end
end

% The perExpt algorithms only identify artifact trials, not time periods.
% Return empty for start
perExpt.start=[];
perExpt.stop = [];
perExpt.trial = find(isArtifactTrial);



end