function flags = detect_bad_behavior_epochs(nsKey,parms)
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
    nsKey (1,1)
    parms.behavior (1,:) string
    parms.bad (1,:) string
    parms.atTrialTime (1,:) double

    parms.exclude = [] % unused
end 

if isa(nsKey, "ns.Experiment")
    expt = nsKey;

else

    expt = ns.Experiment & nsKey;

end

n_bhv = numel(parms.behavior);
flags = ns.NSFlags(parms);
for ii = 1:n_bhv

    statePerTrial  = get(expt,parms.behavior(ii),'prm','state','atTrialTime',parms.atTrialTime(ii),'what','data');
    trials = get(expt,parms.behavior(ii),'prm','state','atTrialTime',parms.atTrialTime(ii),'what','trial');    
    flags.(parms.behavior(ii)) = trials(ismember(upper(statePerTrial),upper(parms.bad(ii))));

end


end