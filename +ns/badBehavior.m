function [perExpt,perChannel] = badBehavior(C,parms)
% Find trials with bad behavior.
% This function is meant to
% be used as the 'fun' in ArtifactParms to fill the Artifact table.
% It is called from populate(ns.Artifact) 
%
% parms.behavior = name of the behavior object (e.g., "fixation")
% parms.bad  = string array identifying the bad behavior(s) (e.g., "FAIL")
% parms.atTrialTime  = time at which to determine the state (e.g., inf for
% end of trial)
%
% EXAMPLE
% fixation.atag = 'fixationBreak';
% fixation.fun = 'ns.badBehavior'; 
% fixation.description  ="Exclude trials with fixation breaks";
% fixation.c = 'eeg';  % Only apply to C files with the eeg ctag
% A bad fixation behavior is defined as the FAIL state at the end of the
% trial
% fixation.parms = struct('behavior','fixation','bad','FAIL','atTrialTime',inf); 
% insertIfNew(ns.ArtifactParm,fixation);

arguments
    C (1,1) ns.C {mustHaveRows(C,1)}
    parms (1,1) struct
end 

expt = ns.Experiment & C; 
statePerTrial  = get(expt,parms.behavior,'prm','state','atTrialTime',parms.atTrialTime,'what','data');
trials = get(expt,parms.behavior,'prm','state','atTrialTime',parms.atTrialTime,'what','trial');

badTrial  = trials(ismember(upper(statePerTrial),upper(parms.bad)));
perExpt.start=[];
perExpt.stop = []; 
perExpt.trial = badTrial(:)';
perChannel = struct([]);  % No per channel artifacts

end