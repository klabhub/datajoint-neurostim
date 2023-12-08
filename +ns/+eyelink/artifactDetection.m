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


%% Get the signals per trial at the sampling rate
[xy,~] = align(C,channel = {'px','py'}, removeArtifacts =false,crossTrial =false,start=parms.start,stop=parms.stop,interpolation="nearest");
[nrSamples,nrTrials,nrChannels] = size(xy); %#ok<ASGLU>

r = sqrt(sum(xy.^2,3));
isArtifactTrial = false(1,nrTrials);
for m = intersect(upper(string(fieldnames(parms))),perExptMethods)'
    switch upper(m)
        case 'SPEED'           
            s = abs(diff(r,1,1));
            outlier = s>parms.speed.max;
            fraction = squeeze(mean(outlier,1,"omitnan")); % Average over samples 
            thisIsArtifactTrial = (fraction>parms.speed.frac);
        case 'ACCELERATION'
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



% Pupil based
[area,~] = align(C,channel = 'pa', removeArtifacts =false,crossTrial =false,start=parms.start,stop=parms.stop,interpolation="nearest");
[nrSamples,nrTrials,nrChannels] = size(area); %#ok<ASGLU>
noChange = diff(area)==0;
fraction = mean(noChange,1);
thisIsArtifactTrial = (fraction>parms.pupil.frac);
fprintf('%d experiment artifact trials based on %s\n',sum(thisIsArtifactTrial),'Pupil');    
isArtifactTrial = isArtifactTrial | thisIsArtifactTrial;

% These algorithms only identify artifact trials, not time periods.
% Return empty for start/
perExpt.start=[];
perExpt.stop = []; 
perExpt.trial = find(isArtifactTrial);


% No per channel artifacts:
perChannel = [] ;



end