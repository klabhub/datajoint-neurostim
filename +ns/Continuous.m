%{
# Preprocessed continuous signals per channel 
-> ns.File
-> ns.ContinuousParm
--- 
startsample: longblob  # Vector of the first sample in each trial 
stopsample: longblob   # Vector with the last sample in each trial
trialstart: longblob   # Start of the trial on the neurostim clock
sampleduration: float  # Duration of a sample
info :  longblob       # General hardware info that is not specific to a channel.
%}
% This is a generic table that (together with its dj.Part table
% PreprocessedChannel) stores preprocessed
% data for electrophysiological recordings. The user provides
% the name of a function that reads data, preprocesses them and returns
% values that are stored in the preprocessed table. Specifically, this
% preprocessing function has the following prototype:
% [signal,time,channelInfo,recordingInfo] = myprep(key, channels, parms)
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
% channelInfo - This is an optional struct array with nrChannels elements. Each
%           struct provides some additional information on the channel to be stored in the
%           info field of the PreprocessedChannel table. For example, this could contain the
%           hardware filtering parameters of the channel as read from the raw data
%           file.
% recordingInfo - General (channel non-specific) info on the recording.
%
% EXAMPLE:
% The ephys.ripple.preprocess functon shows a complete implementation of a
% preprocessing function that handles MUAE, LFP, and EEG recordings with
% the Ripple Grapevine system.
%
% ephys.intan. preprocess does the same for Intan RHX recordings
%
% Once preprocessing is complete, use the get() function of this class to
% retrieve (subsets of ) preprocessed data:
%  For instance, to retrieve the first second in each trial from channels
%  1:10 after preprocessing with the LFP PrepParm:
% [t,v] = get(ephys.Preprocessed & 'prep=''lfp''','channel',1:10,start =0,stop=1)
%
% A sample is assigned to a trial if it occurs after
% the first monitor frame in the trial and before the first monitor frame of the next trial.
% (In other words the ITI is included at the *end* of each trial).
%
% BK - June 2023

classdef Continuous < dj.Computed
    properties (Dependent)
        keySource
    end

    methods

        function v = get.keySource(~)
            % Restrict to files that have the relevant files
            % parms = fetch(ns.ContinuousParm,'extension');
            % restrict = struct('extension',{parms.extension});
            % v =(ns.File & restrict)*proj(ns.ContinuousParm);
            %
            for prms= fetch(ns.ContinuousParm,'extension')'
                restrict  =struct('extension',prms.extension);
                thisV = fetch((ns.File & restrict)*proj(ns.ContinuousParm&restrict));
                if exist("v","var")
                    v  = catstruct(1,v,thisV);
                else
                    v = thisV;
                end
            end
            v = (ns.File*proj(ns.ContinuousParm)) & v;

        end
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
                tbl  (1,1) ns.Continuous
                pv.fetchOptions {mustBeText} = ''
                pv.channel (1,:) double = []
                pv.trial (1,:) double = []
                pv.start (1,1) double = 0
                pv.stop  (1,1) double = 3
                pv.step (1,1) double  = 1/1000;
                pv.interpolation {mustBeText} = 'linear'
                pv.crossTrial (1,1) logical = false;
                pv.align (1,:) double = []
            end

            %% Get the trial mapping
            trialMap = fetch(tbl ,'*');
            if isempty(pv.trial)
                trials = 1:numel(trialMap.trialstart); % All trials
            else
                trials = pv.trial;
            end
            if isempty(pv.align)
                pv.align = zeros(size(trials)); % Align to first frame
            else
                assert(numel(pv.align)==numel(trials),'Each trial should have an align time ')
            end
            % Retrieve the signal in the entire experiment
            % [nrSamples nrTrials]
            if isempty(pv.channel)
                channelRestriction = struct([]);
            else
                channelRestriction = struct('channel',num2cell(pv.channel(:))');
            end
            tblChannel =ns.ContinuousChannel & proj(tbl) & channelRestriction;
            if ~isempty(pv.fetchOptions)
                channelTpl = fetch(tblChannel,'signal',pv.fetchOptions);
            else
                channelTpl = fetch(tblChannel,'signal');
            end
            signal =double([channelTpl.signal]);
            [nrSamples,nrChannels] = size(signal);
            if nrSamples==0||nrChannels==0
                T= timetable; % Empty
            else
                % Setup the new time axis for the results
                newTimes = seconds(pv.start:pv.step:pv.stop)';
                nrTimes  = numel(newTimes);

                % Create a timetable with the activity per trial
                nrTrials = numel(trials);
                varNames = "Trial" + string(trials);
                T =timetable('Size',[nrTimes nrTrials],'RowTimes',newTimes,'VariableTypes',repmat("doublenan",[1 nrTrials]),'VariableNames',varNames);
                trCntr=0;
                % Loop over trials to collect the relevant samples
                trialOut =false(1,nrTrials);
                for tr = trials
                    trCntr= trCntr+1;
                    samplesThisTrial = (trialMap.startsample(tr):trialMap.stopsample(tr))';
                    nrSamplesThisTrial= numel(samplesThisTrial);
                    if isinf(pv.align(trCntr)) || isnan(pv.align(trCntr))
                        fprintf('Align even did not occur in trial %d. Skipping trial.\n',tr)
                        trialOut(trCntr) = true;
                        continue;
                    end
                    trialTime = (0:nrSamplesThisTrial-1)'*trialMap.sampleduration - pv.align(trCntr); %Time in the trial
                    thisT = timetable(seconds(trialTime),signal(samplesThisTrial,:)); % The table for this trial, at the original sampling rate.
                    if pv.crossTrial &&  pv.start <0 && trials(tr) >1
                        % Extract from previous trial (i.e. the time requested was before
                        % firstframe)
                        nrSamplesBefore = ceil(pv.start/trialMap.sampleduration);
                        keepSamples = trialMap.startsample(tr)+(nrSamplesBefore:-1);
                        out = keepSamples <1; % Before the experiment started; remove
                        keepSamples(out) = [];
                        preTime = seconds((-numel(keepSamples):-1)*trialMap.sampleduration); %Time points before current trial start
                        preT = timetable(preTime',signal(keepSamples',:)); % Additional time table with samples before the current trial. (At the original sampling rate)
                        thisT = [preT;thisT]; %#ok<AGROW>
                    end
                    if  pv.crossTrial &&  pv.stop > trialTime(end) && trials(tr) < nrTrials
                        % Extract from next trial (time requested was after
                        % current trial end).
                        nrSamplesAfter = ceil((pv.stop-trialTime(end))/trialMap.sampleduration);
                        keepSamples = trialMap.stopsample(tr)+(1:nrSamplesAfter);
                        out = keepSamples > nrSamples; %After the end of the experiment
                        keepSamples(out) = [];
                        postTime = seconds(trialTime(end) + (1:numel(keepSamples))*trialMap.sampleduration);% Time points after current trial ends
                        postT = timetable(postTime',signal(keepSamples',:)); % Additional time table with samples after the current trial. (At the original sampling rate)
                        thisT = [thisT;postT]; %#ok<AGROW>
                    end
                    % Now retime the table to the new time axis. No
                    % extrapolation
                    thisT = retime(thisT,newTimes,pv.interpolation,'EndValues',NaN);
                    T.(varNames(trCntr)) = table2array(thisT);
                end
                T(:,trialOut) = [];
            end
            % Return as doubles or as timetable.
            if nargout >=2
                [varargout{1},varargout{2}] = timetableToDouble(T);
                if nargout >2
                    nsClockTime = trialMap.trialstart(~trialOut) + repmat(seconds(T.Time),[1 sum(~trialOut)]);
                    varargout{3} = nsClockTime;
                end
            else
                varargout{1} =T;
            end
        end
    end



    methods (Access=protected)
        function makeTuples(tbl,key)
            %% Evaluate the user-specified prep function for the specified parms
            % and channels and insert in the table.
            prepParms = fetch(ns.ContinuousParm & key,'*');            
            % The (user-provided) prep funciton has to return the signal
            % as [nrSamples, nrChannels]  and the channels [nrChannels 1]
            [signal,time,channelInfo,recordingInfo] = feval(prepParms.fun,key,prepParms.parms);
            [nrSamples,nrChannels] = size(signal);
            assert(nrSamples==numel(time),'The number of rows in the preprocessed signal does not match the number of time points ')
            assert(nrChannels==numel(channelInfo),'The number of columns in the preprocessed signal does not match the number of channels')
          
            channels =[channelInfo.nr];

            %% Store the information to map samples to trials (Used in ephys.Preprocessed.get)
            % This is used in the get function to retriev signals per trial
            nrTrials = fetch1(ns.Experiment &key,'trials');
            % Events in neurostim are aligned to firstFrame; we do the same
            % for the analog data
            prms =  get(ns.Experiment & key,'cic');
            trialStartTime = (prms.cic.firstFrameNsTime/1000);
            startSample = nan(1,nrTrials);
            stopSample =nan(1,nrTrials);
            sampleduration = diff(time(1:2));
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
                struct('startsample',startSample, ...
                'stopsample',stopSample,...
                'trialstart',trialStartTime', ...
                'sampleduration',sampleduration, ...
                'info',recordingInfo));
            insert(tbl,tpl)

            % Create tpls for the PreprocessedChannel part table and insert
            channelsTpl = mergestruct(key,...
                struct('signal',num2cell(single(signal),1)',...
                'channel',num2cell(channels(:))));
            if ~isempty(channelInfo)
                for i=1:numel(channelInfo)
                    channelsTpl(i).info = channelInfo(i);
                end
            end

            % Chunking the inserts to avoid overloading the server
            chunkSize = 32; % This should probably be user configurable (e.g., NS_MAXUPLOAD)
            tic;
            fprintf('Uploading to server')
            for i=1:chunkSize:nrChannels
                fprintf('.')
                thisChunk = i:min(nrChannels,i+chunkSize-1);
                insert(ns.ContinuousChannel,channelsTpl(thisChunk));
            end
            fprintf('Done in %d seconds.\n.',round(toc))

        end
    end

    methods (Access=public)
        function [hdr,data,evts] = fieldtrip(tbl,pv)
            arguments
                tbl (1,1) ns.Continuous
                pv.channel =  [];
            end


            %   hdr.Fs                  sampling frequency
            %   hdr.nChans              number of channels
            %   hdr.nSamples            number of samples per trial
            %   hdr.nSamplesPre         number of pre-trigger samples in each trial
            %   hdr.nTrials             number of trials
            %   hdr.label               Nx1 cell-array with the label of each channel
            %   hdr.chantype            Nx1 cell-array with the channel type, see FT_CHANTYPE
            %   hdr.chanunit            Nx1 cell-array with the physical units, see FT_CHANUNIT
            % data - a 2-D matrix of size Nchans*Nsamples for continuous
            % data
            %   event.type      = string
            %   event.sample    = expressed in samples, the first sample of a recording is 1
            %   event.value     = number or string
            %   event.offset    = expressed in samples
            %   event.duration  = expressed in samples
            %   event.timestamp = expressed in timestamp units, which vary over systems (optional)



            for key = fetch(tbl)'
                tpl = fetch(tbl*(proj(ns.ContinuousChannel,'info->channelInfo','signal') &  key) & struct('channel',num2cell(pv.channel)'),'*') ;
                nrChannels = numel(tpl);
                channelInfo = [tpl.channelInfo];
                nrTrials = fetch1(ns.Experiment &key,'trials');
                hdr = struct('Fs',1./tpl.sampleduration,...
                    'nChans',nrChannels, ...
                    'nSamples',size(tpl(1).signal,1), ...
                    'nSamplesPre',0, ...
                    'nTrials',nrTrials,...
                    'label',string(channelInfo.custom_channel_name), ...
                    'chantype',repmat({'eeg'},[1 nrChannels]), ...
                    'chanUnit',repmat({'uv'},[1 nrChannels]));
                %'elec',struct('unit','uv','elecpos','')


                data = [tpl.signal]';

                if nargout>2
                    evts  =struct('type','trial','sample',num2cell(tpl(1).startsample)','value',num2cell(1:numel(tpl(1).startsample))','offset',0,'duration',1,'timestamp',[]);
                end
            end

        end

    end



end