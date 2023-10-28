%{
# Preprocessed continuous signals per channel 
-> ns.File
-> ns.ContParm
--- 
time :      longblob   # Time of each sample (in *milliseconds* on the neurostim clock).
info   :  longblob      # General hardware info that is not specific to a channel.
%}
% This is a generic table that (together with its dj.Part table
% ContChannel) stores data as a continuous vector of samples spanning the
% entire experiment.
%
% The time property stores either all sample times, or in the case of
% regularly sampled data, the time of the first sample, the last sample and
% the number of samples (i.e., [start stop nrSamples]. In the latter case,
% time will be constructed on the fly.
%
% The table is linked to a specific file in an experiment (via ns.File)
% and during a call to populate that file is opened by a function specified
% in the linked ns.ContParm table. This read/preprocess function must have the
% following prototype:
%  [signal,time,channelInfo,recordingInfo] = read(key,parms)
% During a call to populate, key will contain a primary key from the
% ns.File table joined with a key from the ns.ContParm table (i.e. it
% identifies which file to read/preprocess).
% parms will contain the struct stored in ns.ContParm. How this parms
% struct is used is completely up to the read function.(see the example read functions below).
% Before calling the read function, ns.Cont checks that the file exists.
%
% The read function should return the following:
%
% signal -  A matrix with samples along the rows and channels along
%           columns representing all samples obtained during the experiment (ie., all trials).
% time -  The time at which each sample was obtained. For regularly sampled data
%               use a [1 3] vector with [startTime stopTime nrSamples], or irregularly sampled
%               samples, use a [nrSamples 1] column vector. Note that the clock for these times must
%               be in milliseconds and match the neurostim clock. The user function likely
%               includes some code to matchup the clock of the data
%               acquisition device with some event in Neurostim. (see the
%               example read functions below).
% channelInfo - This is a struct array with nrChannels elements. Each
%           struct provides some  information on the channel to be stored in the
%           info field of the ContChannel table. For example, this could contain the
%           hardware filtering parameters of the channel as read from the raw data
%           file. ChannelInfo must contain a .nr field for the number
%           assigned to the channel and can contain a .name field for a
%           descriptive name.
%
% recordingInfo - General (channel non-specific) info on the recording.
%
% EXAMPLE:
% Several read functions have been implemented.
% ephys.intan.read  - Read and Preprocess Intan data
% ns.eyelink.read   - Read EDF files from the Eyelink eye tracker
% ephys.ripple.read - Read and preprocess Ripple Grapevine data.
%
%  To use these, add a row to ContParm with a line like thise
% insert(ns.ContParm,struct('tag','eeg','fun','ephys.intan.read','description'','EEG
%                 data', 'parms',parmsStruct))
% where parmsStruct is a struct that defines what kind of preprocessing the
% read function should do before adding the data to the ns.Cont table.
%
% The align() function of the ns.Cont class is used in the analysis to
% extract trial-based data, or to align to specific events in the
% experiment.
%
% BK - June 2023

classdef Cont< dj.Computed
    properties (Dependent)
        keySource
    end

    methods

        function v = get.keySource(~)
            % Restricted to files with the extenstion specified in ContParm

            % This seems cumbersome, but I coudl not get a simpler join to work
            for prms= fetch(ns.ContParm,'extension','include','exclude')'
                restrict  =struct('extension',prms.extension);
                fileQry = ns.File;
                if ~isempty(prms.include)
                    inc = strsplit(prms.include,',');
                    for i=1:numel(inc)
                          fileQry = fileQry & ['filename LIKE ''' inc{i} ''''];
                    end
                end

                if ~isempty(prms.exclude)
                    exc = strsplit(prms.exclude,',');
                    for i=1:numel(exc)
                          fileQry = fileQry & ['filename NOT LIKE ''' exc{i} ''''];
                    end
                end
                thisV = fetch((fileQry & restrict)*proj(ns.ContParm&restrict));
                if exist("v","var")
                    v  = catstruct(1,v,thisV);
                else
                    v = thisV;
                end
            end
            v = (ns.File*proj(ns.ContParm)) & v;
        end
    end
    methods (Access=public)
        function [varargout] = align(tbl,pv)
            % Function to retrieve trial-based and aligned preprocessed signals.
            %
            % A sample is assigned to a trial if it occurs after the first monitor frame
            % in the trial and before the first monitor frame of the next trial.
            % (In other words the ITI is included at the *end* of each trial).
            % INPUT
            % tbl  - ns.Cont table
            %
            % Optional Parameter/Value pairs
            % channel  - The subset of channels ([]; means all)
            % trial   - Which trials to extract
            % The temporal snippet of data to extract runs from start to
            % stop in steps of step.These times are specified relative to
            % the align event (below) in milliseconds.
            % start  - Extract signals startint at this time (ms)
            % stop   - Last time point to extract. (ms)
            % step   - Step size in seconds. (ms)
            % align - The time in each trial that is considered 0. By
            % default this is the time at which the first monitor frame
            % became visible. By specifying the time of an event (e.g.,
            % stimulus onset), the data can be aligned to that event.
            %
            % interpolation -  enum('nearest','linear','spline','pchip','makima')
            %               Interpolation method; see timetable/synchronize. ['linear']
            % crossTrial - Allow values to be returned that are from the
            %               trials before or after. When this is false,
            %               only times/samples between the firstFrame event in the trial and
            %               the firstFrame event of the next trial will be
            %               returned.
            % fetchOptions - Options passed to the fectch(ns.Cont)
            %               call.(For instance 'LIMIT 1')
            %
            % OUTPUT
            %  [v,t]  = t: time in milliseconds since first frame event,
            %           v: Matrix with [nrTimePoints nrTrials nrChannels]
            % Alternatively, when only a single output is requested:
            % T     = timetable with each column a trial. Time is in seconds
            %           relative to the first frame of the trial.
            %          Channels are along the columns of the rows of the
            %           elements of the table. The Time dimension will be 
            %           named TrialTime  and the data dimension named after
            %               the preprocessing tag.
            %
            % EXAMPLE
            % T= align(ns.Cont & 'tag=''eeg''',channel=2,align=
            % stimulusOnsetTime, start =-250,stop=1000,step=1)
            % will return a timetable with the EEG samples of channel 2,
            % from 250 ms before until 1000 ms after the stimulus onset
            % time in each trial. (stimulusOnsetTime has to be a vector
            % with a length equal to the number of trials in which each
            % element represent the onset time in ms such as returned by ns.Experiment.get()
            arguments
                tbl  (1,1) ns.Cont {mustHaveRows(tbl,1)}
                pv.fetchOptions {mustBeText} = ''
                pv.channel (1,:) double = []
                pv.trial (1,:) double = []
                pv.start (1,1) double = 0
                pv.stop  (1,1) double = 1000
                pv.step (1,1) double  = 1;
                pv.interpolation {mustBeText} = 'linear'
                pv.crossTrial (1,1) logical = false;
                pv.align (1,:) double = []
            end

            %% Collect data from dbase and check consistency
            sampleTime = fetch1(tbl ,'time');
            if numel(sampleTime)==3
                sampleTime = linspace(sampleTime(1),sampleTime(2),sampleTime(3))';
            end
            nrAllTrials = fetch1(ns.Experiment*tbl,'trials');
            allTrials = 1:nrAllTrials;
            if isempty(pv.trial)
                trials = allTrials; % All trials
            else
                trials = intersect(pv.trial,allTrials);
            end

            if isempty(pv.align)
                pv.align = zeros(size(trials)); % Align to first frame
            elseif numel(pv.align)==1 % Singleton expansion
                pv.align = repmat(pv.align,[1 numel(trials)]);
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
            tblChannel =ns.ContChannel & proj(tbl) & channelRestriction;
            if ~isempty(pv.fetchOptions)
                channelTpl = fetch(tblChannel,'signal',pv.fetchOptions);
            else
                channelTpl = fetch(tblChannel,'signal');
            end
            signal =double([channelTpl.signal]); % Signal as matrix
            [nrSamples,nrChannels] = size(signal);
            % Fetch the preprocessing tag. The dimensions in the timetable
            % will be named after this.
            tag = fetch1(ns.ContParm & tbl,'tag');
            %% Align
            if nrSamples==0||nrChannels==0
                T= timetable; % Empty
            else
                % Setup the new time axis for the results
                newTimes = milliseconds(pv.start:pv.step:pv.stop)';
                nrTimes  = numel(newTimes);

                % Create a timetable with the activity per trial
                nrTrials = numel(trials);
                varNames = "Trial" + string(trials);
                T =timetable('Size',[nrTimes nrTrials],'RowTimes',newTimes,'VariableTypes',repmat("doublenan",[1 nrTrials]),'VariableNames',varNames);
                % trial time zero in NS for each trial
                trialStartTime = get(ns.Experiment & tbl, 'cic','prm','firstFrame','atTrialTime',inf,'what','clocktime');

                % Loop over trials to collect the relevant samples
                trialOut =false(1,nrTrials);

                for trCntr=1:nrTrials
                    if isinf(pv.align(trCntr)) || isnan(pv.align(trCntr))
                        fprintf('Align even did not occur in trial %d. Skipping trial.\n',trials(trCntr))
                        trialOut(trCntr) = true;
                        continue;
                    end
                    thisTrial =trials(trCntr);
                    % Limits if crossTrial is allowed
                    startLimit =  trialStartTime(thisTrial) + pv.align(trCntr)+pv.start;
                    stopLimit =  trialStartTime(thisTrial) + pv.align(trCntr)+pv.stop; 
                    if ~pv.crossTrial
                        startLimit =  max(trialStartTime(thisTrial),startLimit); % No samples before trialStartTime(thisTrial)
                        if thisTrial<nrAllTrials
                            stopLimit = min(trialStartTime(thisTrial+1),stopLimit); % No sample after the start of next trial (last trial includes everything until the end of recording).
                        end
                    end
                    staySamples = sampleTime >= startLimit & sampleTime < stopLimit;
                    
                    % Unless the first stay sample is exactly at
                    % trialstart, the first sample in the new table will be
                    % NaN. To circumvent this, add one sample before (and
                    % after) the range we identified.
                    ixFirst = find(staySamples,1,'first');                   
                    if ixFirst>1
                        staySamples(ixFirst-1) =true;
                    end
                    ixLast= find(staySamples,1,'last');                   
                    if ixLast<nrSamples
                        staySamples(ixLast+1) =true;
                    end

                    trialTime = sampleTime(staySamples)-trialStartTime(thisTrial)- pv.align(trCntr); %Time in the trial
                    thisT = timetable(milliseconds(trialTime),signal(staySamples,:)); % The table for this trial, at the original sampling rate.
                    % Now retime the table to the new time axis. Never extrapolation
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
                T.Properties.DimensionNNames = {'TrialTime',tag};
            end
        end
    end



    methods (Access=protected)
        function makeTuples(tbl,key)

            % Check that the relevant file exists
            qry = ns.File & key;
            nrFiles = count(qry);
            if nrFiles ~=1
                % Zero or more than 1 file
                error('This experiment has %d files. Cannot proceed.',nrFiles);
            else
                % Fetch the file to read
                filename = fullfile(folder(ns.Experiment &key),fetch1(qry,'filename'));
            end
            if exist(filename,'file') && ~exist(filename,'dir')
                fprintf('Reading %s\n',filename);
            else
                error('File %s does not exist. Cannot create Cont',filename);
            end

            % Get the specified preprocessing parms
            prepParms = fetch(ns.ContParm & key,'*');
            % Call the prep function
            [signal,time,channelInfo,recordingInfo] = feval(prepParms.fun,key,prepParms.parms);
            [nrSamples,nrChannels] = size(signal);

            assert(nrSamples==numel(time) || (numel(time)==3 &&time(3)==nrSamples),'The number of rows in the preprocessed signal does not match the number of time points ')
            assert(nrChannels==numel(channelInfo),'The number of columns in the preprocessed signal does not match the number of channels')

            channels =[channelInfo.nr];

            % Create tuples and insert.
            tpl = mergestruct(key,...
                struct('time',time, ...
                'info',recordingInfo));
            insert(tbl,tpl)

            % Create tpls for each of the channels and insert
            channelsTpl = mergestruct(key,...
                struct('signal',num2cell(single(signal),1)',...
                'channel',num2cell(channels(:))));
            if ~isempty(channelInfo)
                for i=1:numel(channelInfo)
                    channelsTpl(i).info = channelInfo(i);
                    if isfield(channelInfo,'name')
                        channelsTpl(i).name = channelInfo(i).name;
                    end
                end
            end

            % Chunking the inserts to avoid overloading the server
            chunkSize = 32; % This should probably be user configurable (e.g., NS_MAXUPLOAD)
            tic;
            fprintf('Uploading to server ')
            for i=1:chunkSize:nrChannels
                fprintf('.')
                thisChunk = i:min(nrChannels,i+chunkSize-1);
                insert(ns.ContChannel,channelsTpl(thisChunk));
            end
            fprintf('Done in %d seconds.\n.',round(toc))

        end
    end

    methods (Access=public)
        function [hdr,data,evts] = fieldtrip(tbl,pv)
            arguments
                tbl (1,1) ns.Cont
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
                tpl = fetch(tbl*(proj(ns.ContChannel,'info->channelInfo','signal') &  key) & struct('channel',num2cell(pv.channel)'),'*') ;
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