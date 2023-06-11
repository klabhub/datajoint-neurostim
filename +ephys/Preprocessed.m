%{
# Preprocessed signals per electrode and experiment 
-> ns.Experiment
-> ephys.PrepParm
--- 
%}
% This is a generic table that (together with its dj.Part tables
% PreprocessedChannel and PreprocessedTrialmap stores preprocessed
% data for electrophysiological recordings. The user provides
% the name of a function that reads data, preprocesses them and returns
% values that are stored in the preprocessed table. Specifically, this
% preprocessing function has the following prototype:
% [signal,time,info] = myprep(key, channels, parms)
% INPUT
% key   -  The tuple that defines the source data (i.e. a row in the join of ns.Experiment
%               and ephys.PrepParm). The user function uses this key
%               to load the appropriate file  (e.g. from ns.File & key),
%
% channels -  The channels to preprocess. This input to the function is
%               determined from the Array field of the PrepParms table
%               (ephys.Array & key).
% parms  - A struct with parameters (from (ephys.PrepParm &key)). The user
%           function uses this to determine what kind of preprocessing to do.
% OUTPUT
% signal -  A matrix with samples along the rows and channels along
%           columns. These columns must match the channels input to the function,
%           This contains all samples obtained during the experiment (ie., all trials).
% time -  The time at which each sample was obtained. This is a [nrSamples 1]
%               column vector. Note that the clock for these times must
%               match the neurostim clock. The user function likely
%               includes some code to matchup the clock of the data
%               acquisition device with some event in Neurostim.
% info - This is an optional struct array with nrChannels elements. Each cell
%           provides some additional information on the channel that is stored in the
%           info field of the Preprocessed table. For example, this could contain the
%           hardware filtering parameters of the channel as read from the raw data
%           file.
%
% EXAMPLE:
% The ephys.ripple.preprocess functon shows a complete implementation of a
% preprocessing function that handles MUAE, LFP, and EEG recordings with
% the Ripple Grapevine system.
% 
% Once preprocessing is complete, use the get() function of this class to
% retrieve (subsets of ) preprocessed data:
%  For instance, to retrieve the first second in each trial from channels
%  1:10 after preprocessing with the LFP PrepParm:
% [t,v] = get(ephys.Preprocessed & 'prep=''lfp''','channel',1:10,start =0,stop=1)
%
% BK - June 2023

classdef Preprocessed < dj.Imported
    properties (Dependent)
        %    keySource  (limit to a subset of experiments with ephys data)
    end

    methods (Access=public)
        function [varargout] = get(tbl,pv)
            % Function to retrieve trial-start aligned preprocessed signals
            % per channel.
            %
            % tbl  - ephys.Preprocessed table /rows to use
            %
            % Optional Parameter/Value pairs
            % trial   - Which trials to extract
            % start  - First time point to extract (relatve to first frame
            %           of each trial, in seconds)
            % stop   - Last time point to extract.
            % step   - Step size in seconds.
            % interpolation -  enum('nearest','linear','spline','pchip','makima')
            %               Interpolation method; see timetable/synchronize. ['linear']
            % crossTrial - Allow values to be returned that are from the
            %               trials before or after. When this is false,
            %               only times/samples between the firstFrame event in the trial and
            %               the firstFrame event of the next trial will be
            %               returned.
            %           
            % OUTPUT
            %  [v,t]  = t: time in seconds since first frame event,
            %           v: Matrix with [nrTimePoints nrTrials nrChannels]
            % Alternatively, when only a single output is requested:
            % T     = timetable with each column a trial. Time is in seconds
            %           relative to the first frame of the trial.
            %          Channels are along the columns of the rows of the 
            %           elements of the table.           
            arguments
                tbl  (1,1) ephys.Preprocessed {mustHaveRows}
                pv.fetchOptions {mustBeText} = ''
                pv.channel (1,:) double = []
                pv.trial (1,:) double = []
                pv.start (1,1) double = 0
                pv.stop  (1,1) double = 3
                pv.step (1,1) double  = 1/1000;
                pv.interpolation {mustBeText} = 'linear'
                pv.crossTrial (1,1) logical = false;
            end

            %% Get the trial mapping
            trialMap = ephys.PreprocessedTrialmap & tbl;
            trialMap = fetch(trialMap ,'*');
            if isempty(pv.trial)
                trials = [trialMap.trial]; % All trials
            else
                trials = pv.trial;
            end
            % Retrieve the signal in the entire experiment
            % [nrSamples nrTrials]
            if isempty(pv.channel)
                channelRestriction = struct([]);
            else
                channelRestriction = struct('channel',num2cell(pv.channel(:))');
            end
            tblChannel =ephys.PreprocessedChannel & tbl & channelRestriction;
            if ~isempty(pv.fetchOptions)
                channelTpl = fetch(tblChannel,'signal',pv.fetchOptions);
            else
                channelTpl = fetch(tblChannel,'signal');
            end
            signal =[channelTpl.signal];
            [nrSamples,nrChannels] = size(signal); %#ok<ASGLU>
            % Setup the new time axis for the results
            newTimes = seconds(pv.start:pv.step:pv.stop)';
            nrTimes  = numel(newTimes);

            % Create a timetable with the activity per trial
            sampleDuration = unique([trialMap.sampleduration]);
            nrTrials = numel(trials);
            varNames = "Trial" + string([trialMap(trials).trial]);
            T =timetable('Size',[nrTimes nrTrials],'RowTimes',newTimes,'VariableTypes',repmat("doublenan",[1 nrTrials]),'VariableNames',varNames);
            trCntr=0;
            % Loop over trials to collect the relevant samples
            for tr = trials
                trCntr= trCntr+1;
                samplesThisTrial = (trialMap(tr).startsample:trialMap(tr).stopsample)';
                nrSamplesThisTrial= numel(samplesThisTrial);
                trialTime = (0:nrSamplesThisTrial-1)'*sampleDuration; %Time in the trial
                thisT = timetable(seconds(trialTime),signal(samplesThisTrial,:)); % The table for this trial, at the original sampling rate.
                if pv.crossTrial &&  pv.start <0 && trialMap(tr).trial>1
                    % Extract from previous trial (i.e. the time requested was before
                    % firstframe)
                    nrSamplesBefore = ceil(pv.start/sampleDuration);
                    keepSamples = trialMap(tr).startsample+(nrSamplesBefore:-1);
                    out = keepSamples <1; % Before the experiment started; remove
                    keepSamples(out) = [];
                    preTime = seconds((-numel(keepSamples):-1)*sampleDuration); %Time points before current trial start
                    preT = timetable(preTime',signal(keepSamples',:)); % Additional time table with samples before the current trial. (At the original sampling rate)
                    thisT = [preT;thisT]; %#ok<AGROW>
                end
                if  pv.crossTrial &&  pv.stop > trialTime(end) && trialMap(tr).trial < numel(trialMap)
                    % Extract from next trial (time requested was after
                    % current trial end).
                    nrSamplesAfter = ceil((pv.stop-trialTime(end))/sampleDuration);
                    keepSamples = trialMap(tr).stopsample+(1:nrSamplesAfter);
                    out = keepSamples > nrSamples; %After the end of the experiment
                    keepSamples(out) = [];
                    postTime = seconds(trialTime(end) + (1:numel(keepSamples))*sampleDuration);% Time points after current trial ends
                    postT = timetable(postTime',signal(keepSamples',:)); % Additional time table with samples after the current trial. (At the original sampling rate)
                    thisT = [thisT;postT]; %#ok<AGROW>
                end
                % Now retime the table to the new time axis. No
                % extrapolation
                thisT = retime(thisT,newTimes,pv.interpolation,'EndValues',NaN);
                T.(varNames(trCntr)) = table2array(thisT);
            end
            % Return as doubles or as timetable.
            if nargout ==2
                [varargout{1},varargout{2}] = timetableToDouble(T);
            else
                varargout{1} =T;
            end
        end
    end


    methods (Access=protected)
        function makeTuples(tbl,key)

            %% Inser the  key in the Preprocessed table. This is just a
            % placeholder.
            insert(tbl,key)

            %% Evaluate the user-specified prep function for the specified parms
            % and channels and insert in the table.
            preParms = fetch(ephys.PrepParm & key,'*');
            channels = fetch1(ephys.Array & key,'channels');
            % The (user-provided) prep funciton has to return the signal
            % as [nrSamples, nrChannels]  and the channels [nrChannels 1]
            [signal,time,info] = feval(preParms.fun,key,channels, preParms.parms);
            [nrSamples,nrChannels] = size(signal);
            assert(nrSamples==numel(time),'The number of rows in the preprocessed signal does not match the number of time points ')
            assert(nrChannels==numel(channels),'The number of columns in the preprocessed signal does not match the number of channels')
            assert(isempty(info) || (numel(info)==nrChannels),'The info rerturned by preprocessing deos not match the number of channels')

            % Create tpls for the PreprocessedChannel part table and insert
            channelsTpl = mergestruct(key,...
                struct('signal',num2cell(signal,1)',...
                'channel',num2cell(channels(:))));
            if ~isempty(info)
                for i=1:numel(info)
                    channelsTpl(i).info = info(i);
                end
            end
            insert(ephys.PreprocessedChannel,channelsTpl);


            %% Setup a table to map samples to trials (Used in ephys.Preprocessed.get)
            % This is used in the get function to retriev signals per trial
            nrTrials = fetch1(ns.Experiment &key,'trials');
            % Events in neurostim are aligned to firstFrame; we do the same
            % for the analog data
            prms =  get(ns.Experiment & key,'cic');
            trialStartTime = (prms.cic.firstFrameNsTime/1000);
            startSample = nan(1,nrTrials);
            stopSample =nan(1,nrTrials);
            sampleDuration = diff(time(1:2));
            for tr=1:nrTrials
                stay = time >=trialStartTime(tr);
                if tr <nrTrials
                    stay = stay & time <trialStartTime(tr+1);
                end
                startSample(tr) = find(stay,1,'first');
                stopSample(tr) = find(stay,1,'last');
            end
            % Create tuples and insert.
            tpl = mergestruct(key,...
                struct('trial',num2cell(1:nrTrials)',...
                'startsample',num2cell(startSample)', ...
                'stopsample',num2cell(stopSample)',...
                'trialstart',num2cell(trialStartTime), ...
                'sampleduration',repmat({sampleDuration},[nrTrials 1])));

            insert(ephys.PreprocessedTrialmap,tpl);

        end
    end

end