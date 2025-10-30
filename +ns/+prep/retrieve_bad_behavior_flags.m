function flags = retrieve_bad_behavior_flags(nsKey,parms)
% ns.prep.retrieve_bad_behavior_flags
% previously known as ns.detect_bad_behavior_epochs
% Find trials with bad behavior, flagged during recording
% This function is meant to
% be used as the 'fun' in artifact_parm of EpochParm
%
% parms.behavior = name of the behavior object (e.g., "fixation")
% parms.bad  = string array identifying the bad behavior(s) (e.g., "FAIL")
% parms.atTrialTime  = time at which to determine the state (e.g., inf for
% end of trial)
%
% EXAMPLE
% args = {'behavior','fixation','bad','FAIL','atTrialTime',inf};
% flags = ns.prep.retrieve_bad_behavior_flags(key, args{:})
% returns flags for fixation breaks

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
flags = ns.prep.NSFlags(parms);
for ii = 1:n_bhv

    statePerTrial  = get(expt,parms.behavior(ii),'prm','state','atTrialTime',parms.atTrialTime(ii),'what','data');
    trials = get(expt,parms.behavior(ii),'prm','state','atTrialTime',parms.atTrialTime(ii),'what','trial');    
    flags.(parms.behavior(ii)) = trials(ismember(upper(statePerTrial),upper(parms.bad(ii))));

end


end