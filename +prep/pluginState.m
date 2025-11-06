function BB = pluginState(expt,trials,parms)
% Select trials on the basis of plugin properties. This is used by
% ns.Epoch to remove trials that failed (e.g. using fixation events)
%  Returns a prep.badBy object with the bad trials.
%
% The parms struct defines the good trials.
%  EXAMPLE
% parms.fix.plugin = "fixation"         -> look in the plugin named fixation
% parms.fix.parameter = "state"         -> consider the state parameter
% parms.fix.atTrialTime = Inf            -> at the end of the trial
% parms.fix.mode = "data"               -> read the data
% parms.fix.data = "SUCCESS"            -> and compare those data with the  value "SUCCESS"
%
% Additional criteria can be added by adding more fields to the parms
% struct. (The name of the field is only used for reporting).
%
% This will select trials in which the 'state' parameter  of the fixation
% plugin  has the value SUCCESS at the end of the trial.
arguments
    expt (1,1) ns.Experiment
    trials (1,:) double
    parms (1,1) struct {prep.mustBePlgParm}  % Struct with instructtions
end
BB = prep.badBy;
fn = fieldnames(parms);
nrCriteria= numel(fn);
for c= 1:nrCriteria  % Loop over all criteria
    BB.(fn{c}) = []; % Initialize empty
    thisParms = parms.(fn{c});
    if ~isfield(thisParms,'mode')
        thisParms.mode ="data"; % Default to comparing data
    end
    switch thisParms.mode
        case "data"
            plgData = get(expt,thisParms.plugin,'prm',thisParms.parameter,'atTrialTime',thisParms.atTrialTime,'what',["data" "trial"],'trial',trials);
            if ischar(thisParms.data) || isstring(thisParms.data) || iscellstr(thisParms.data)
                ok = ismember(string(plgData.data),string(thisParms.data));
            else
                ok = ismember(plgData.data,thisParms.data);
            end
            BB.(fn{c}) = plgData.trial(~ok);
            % case "time" - could add to check whether an event occurs at a
            % certain time.
        otherwise
            error("niy")
    end
end
end