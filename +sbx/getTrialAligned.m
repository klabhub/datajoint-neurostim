function [varargout] = getTrialAligned(tbl,expt,data,pv)
% Function to retrieve trial-start aligned data for a Eye or Ball table.
% tbl  - Table to use
% key - Specific experiment (tuple from fetch)
% data - name of the data column to return (e.g. 'x','y',z' in sbx.Eye)
%
% Optional Parameter/Value pairs
% trial   - Which trials to extract
% start  - First time point to extract (relatve to first frame
%           of each trial, in seconds)
% stop   - Last time point to extract
% step   - Step size in seconds.
% interpolation -  enum('nearest','linear','spline','pchip','makima')
%               Interpolation method; see timetable/retime. ['linear']
%                   For binning, use one of the aggregation methods
%                   ('mean','median')
% crossTrial - Allow values to be returned that are from one
% trial before or one trial after. This is helpful to set start
% =-1 to get the values from the iti before the trial. [true]
%
% OUTPUT
% T     = timetable with each column a trial. Time is in seconds
%           relative to the first frame of the trial. The different data
%           columns are columns per trial.
% To return one timetable per data column, specify a matching number of
% output arguments.
%
arguments
    tbl (1,1)
    expt (1,1) ns.Experiment {mustHaveRows(expt,1)}
    data (1,:) cell

    pv.trial (1,:) double = []
    pv.start (1,1) double = 0
    pv.stop  (1,1) double = 3
    pv.step (1,1) double  = 0.250;
    pv.interpolation {mustBeText} = 'linear'
    pv.crossTrial (1,1) logical = true;
end

%% Get the mapping from Frames to trials.
trialMap = fetch(ns.MovieTrialmap & tbl & expt,'*');
noFrames= cellfun(@ischar,{trialMap.frame});
trialMap(noFrames) = [];
if isempty(pv.trial)
    trials = [trialMap.trial]; % All trials
else
    trials = pv.trial;
end
%% Fetch the data for the full experiment.
nrDataColumns = numel(data);
if nargout>nrDataColumns
    error('Number of outputs (%d) should not be larger than the number of requested data columns (%d)',nargout,nrDataColumns);
end
V = cell(1,nrDataColumns);
[V{:}]= fetchn(tbl&expt,data{:}); % fetchn wraps in a cell
V = cellfun(@(x)(x(1)),V);
V = [V{:}];  % One column per column...

%% Setup the target table
newTimes = seconds(pv.start:pv.step:pv.stop);
nrTimes  = numel(newTimes);
frameDuration = seconds(1./unique([fetch(ns.Movie*(tbl& expt),'framerate').framerate]));
nrTrials = numel(trials);
varNames = "Trial" + string([trialMap(trials).trial]);
T =timetable('Size',[nrTimes nrTrials],'RowTimes',newTimes','VariableTypes',repmat("doublenan",[1 nrTrials]),'VariableNames',varNames);
trCntr=0;

%% Loop over the trials to extract the relevant frames from the movie-derived data.
for tr = trials
    trCntr= trCntr+1;
    thisT = timetable(seconds(trialMap(tr).trialtime(:)),V(trialMap(tr).frame,:));

    % With Crosstrial set to true we take frames from preceding or
    % succeding trials, assuming that the frameDuration is constant (i.e.
    % regular sampling by the TPI). 
    if pv.crossTrial &&  pv.start <0 && trialMap(tr).trial>1
        % Extract from previous trial (i.e. the time requested was before
        % firstframe)
        nrFramesBefore = ceil(pv.start/seconds(frameDuration));
        framesBefore  =nrFramesBefore:-1;
        preTime  = seconds(trialMap(tr).trialtime(1))+framesBefore*frameDuration;
        keepFrames = trialMap(tr).frame(1)+framesBefore;
        out = keepFrames <1;
        preTime(out) = [];
        keepFrames(out) = [];
        preT = timetable(preTime',V(keepFrames',:));
        thisT = [preT;thisT]; %#ok<AGROW>
    end

    if  pv.crossTrial &&  pv.stop > trialMap(tr).trialtime(end) && trialMap(tr).trial < numel(trialMap)
        % Extract from next trial (time requested was after trial end).
        nrFramesAfter = ceil((pv.stop-trialMap(tr).trialtime(end))/seconds(frameDuration));
        framesAfter =  (1:nrFramesAfter);
        postTime     = seconds(trialMap(tr).trialtime(end))+framesAfter*frameDuration;
        keepFrames = trialMap(tr).frame(end)+framesAfter;
        out = keepFrames >trialMap(end).frame(end);
        postTime(out) = [];
        keepFrames(out) = [];
        postT = timetable(postTime',V(keepFrames',:));
        thisT = [thisT;postT]; %#ok<AGROW>
    end
    thisT = retime(thisT,newTimes,pv.interpolation,'EndValues',NaN);
    T.(varNames(trCntr)) = table2array(thisT);
end
aggregationMethods= {'sum','mean','median','mode','prod','min','max','firstValue','lastValue'};
if ismember(pv.interpolation,aggregationMethods)
    % retime with aggregation returns the left edge of each bin
    % and the last entry in the table is the value that occurs
    % exactly at the last time point (so not a bin).  Correct
    % this here to return the time of the center of the bins
    % and only the mean values in those bins.
    T(end,:) = [];
    T.Time = T.Time -seconds(pv.step)/2;
end

%% Prepare the output
if nargout ==1
    % Single timetable
    varargout{1} =T;
elseif nargout == nrDataColumns
    % One timetable per data column
    varargout=  cell(nargout,1);
    vals= T.Variables;
    for i=1:nrDataColumns
        varargout{i} = array2timetable(vals(:,i:nrDataColumns:end),'RowTimes',T.Time,'VariableNames',varNames);
    end
end
end