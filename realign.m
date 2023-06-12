function T = realign(T,newTime)
% Realign a timetable to a new set of times
% INPUT
% T = A timetable as returned by ephys.Preprocessed get. It has samples
% along rows, and trials along the columns. 
% newTime = The new time to align to. when this is a scalar it merely
% redefines the T.Time, when this is a vector with one entry for each
% trial, each trial is realigned to that event.
% OUTPUT
% T - The timetable with each trial aligned to the new time.
% 
% EXAMPLE:
% if an event (say stimulus onset) happend 1 s into the first trial and 2 s
% into the second trial, then we can align the data in a T with two trials 
% to stimulus onset with T = realign(T,[1 2]);
arguments 
    T (:,:) timetable
    newTime (1,:) double
end

nrTimes = numel(newTime);
assert(nrTimes==1 || nrTimes==width(T),'The number of times for realign should be a singleton or 1 per trial');
[~,nrTrials] =size(T);
if nrTimes ==1
    % Only redefinig the time axis; same for all trials
    % This is called recursively for each trial.
    % Find the nearest t
    nrTimeSteps = round(seconds(newTime)/T.Properties.TimeStep);
    shift = nrTimeSteps*T.Properties.TimeStep;
    T.Time = T.Time -shift;
else
    % Each trial has its own time.
    newT = realign(T(:,1),newTime(1));
    for tr =2:nrTrials
        newT = synchronize(newT,realign(T(:,tr),newTime(tr)));
    end
    T= newT;        
end