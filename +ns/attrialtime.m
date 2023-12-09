function out = attrialtime(props,propName,time,c,what,trial)
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
% what - 'data','trialtime','clocktime','trial': which aspect of the data
%               to return.
% trial - a vector of trial numbers for which to return the information.
% Defaults to [], which means all trials.
%
% OUTPUT
% values  = A cell array of values, one for each requested trial.
%
% BK - Dec 2022.
if nargin <6
    trial = [];
    if nargin <5
        what = 'data';
    end
end

if isnan(time)
    switch what
        case 'data'
            out = props.(propName);
        case 'trialtime'
            out = props.([propName 'Time']);
        case 'clocktime'
            out = props.([propName 'NsTime']);
        case 'trial'
            out = props.([propName 'Trial']);
    end
    if ~isempty(trial)
        trialNrs = props.([propName 'Trial']);
        out = out(ismember(trialNrs,trial));
    end
    return;
end

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

if isfield(props,[propName 'NsTime'])
    allEventNsTimes  = props.([propName 'NsTime']);
else
    allEventNsTimes = -inf;
end




% Initialize with nan
data = cell(nrTrials,1);
eventTime = nan(nrTrials,1);
eventNsTime = nan(nrTrials,1);
[data{:}] = deal(NaN);
currentTrial =NaN;
% Loop through all events (= events when the property changed value)
% For many events this could be made more efficient by skipping successive
% events in a given trial, but the find needed for this is probably slower
% than this loop over all events
for e=1:numel(allEventTrials)
    if ~isnan(currentTrial) && allEventTrials(e) > currentTrial+1
        % Fill in skipped trials with the currentValue
        trgTrials = currentTrial+1:allEventTrials(e)-1;
        [data{trgTrials}] =deal(currentValue);
        eventTime(trgTrials) = -inf; % Indicates that this value was set before the start of the trial
        eventNsTime(trgTrials) = currentNsTime;
    end
    % Next trial or same trial, update until atTrialTime reached.
    currentTrial = allEventTrials(e);
    currentTime = allEventTimes(e);
    currentNsTime = allEventNsTimes(e);
    if iscell(allEventValues)
        currentValue = allEventValues{e};
    else
        currentValue = allEventValues(e);
    end
    if currentTime <= time
        data{currentTrial} = currentValue;
        eventTime (currentTrial) = currentTime;
        eventNsTime(currentTrial) = currentNsTime;
    else
        % Event occurred after the atTrialTime; use the value from the
        % previous trial
        if currentTrial >1
            data{currentTrial} = data{currentTrial-1};
            eventTime (currentTrial) = -inf;
            eventNsTime(currentTrial) = eventNsTime(currentTrial-1);
        end
    end
end
% Fill in to the end from the last value that was stored.
if currentTrial~=nrTrials
    [data{currentTrial+1:nrTrials}] =deal(currentValue);
    eventNsTime(currentTrial+1:nrTrials) = currentNsTime;
    eventTime(currentTrial+1:nrTrials)  =-inf;
end

switch what
    case 'data'
        out = data ;
    case 'trialtime'
        out = eventTime;
    case 'clocktime'
        out = eventNsTime;
    case 'trial'
        out = find(~isinf(eventNsTime)); % Trials in whcih the event actually ocurred
end

if ~isempty(trial)
        trialNrs = find(~isinf(eventNsTime));
        out = out(ismember(trialNrs,trial));
end
    
end