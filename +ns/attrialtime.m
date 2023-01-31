function vals = attrialtime(props,propName,time,c)
% Given a struct with properties (as returned by ns.Experiment.get  or
% nsPluginParameter.get, and a specific property, return the value of that
% property at a certain time in each trial. 
% This function is necessary to deal with the "feature" that Neurostim
% stores only changes to a property, hence if a property did not change in
% trial n, it shoudl have the value it got in trial n-1. 
%
% INPUT
% props - struct with properties , as returned by one ofthe get functions.
% propName - The name of the property that shoudl be returned.
% time - The time (relative to the trial start) at which the property
%           should be determined. Use Inf for 'at the end of the trial'.
% c - A struct with the properties of the CIC (also returned by
%               ns.Experiment.get)
% OUTPUT
% values  = A cell array of values, one for each trial. 
% 
% BK - Dec 2022.


nrTrials  = max(c.trial);
allEventValues = props.(propName);
if isfield(props,[propName 'Trial'])
allEventTrials = props.([propName 'Trial']);
else
    allEventTrials =1;
end
if isfield(props,[propName 'Time'])
    allEventTimes  = props.([propName 'Time']);
else
    allEventTimes = -inf;
end

% Initialize with nan
vals = cell(nrTrials,1);
[vals{:}] = deal(NaN);
currentTrial =NaN;
% Loop through all events (= events when the property changed value)
% For many events this could be made more efficient by skipping successive
% events in a given trial, but the find needed for this is probably slower
% than this loop over all events
for e=1:numel(allEventTrials)
    if ~isnan(currentTrial) && allEventTrials(e) > currentTrial+1
        % Fill in skipped trials with the currentValue
        [vals{currentTrial+1:allEventTrials(e)-1}] =deal(currentValue);
    end    
    % Next trial or same trial, update until atTrialTime reached.
    currentTrial = allEventTrials(e);
    currentTime = allEventTimes(e);
    if iscell(allEventValues)
        currentValue = allEventValues{e};
    else
        currentValue = allEventValues(e);
    end
    if currentTime <= time
        vals{currentTrial} = currentValue;
    end    
end
% Fill in to the end from the last value that was stored.
if currentTrial~=nrTrials
     [vals{currentTrial+1:nrTrials}] =deal(currentValue);
end
