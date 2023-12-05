%{
# Preprocessed continuous signals per channel 
-> ns.File
-> ns.CParm
--- 
time :      longblob   # Time of each sample (in *milliseconds* on the neurostim clock).
info   :  longblob      # General hardware info that is not specific to a channel.
nrsamples: int          # Number of samples/frames across the experiment
%}
% This is a generic table that (together with its dj.Part table
% CChannel) stores data as a continuous vector of samples spanning the
% entire experiment.
%
% The time property stores either all sample times, or in the case of
% regularly sampled data, the time of the first sample, the last sample and
% the number of samples (i.e., [start stop nrSamples]. In the latter case,
% time will be constructed on the fly.
%
% The table is linked to a specific file in an experiment (via ns.File)
% and during a call to populate that file is opened by a function specified
% in the linked ns.CParm table. This read/preprocess function must have the
% following prototype:
%  [signal,time,channelInfo,recordingInfo] = read(key,parms)
% During a call to populate, key will contain a primary key from the
% ns.File table joined with a key from the ns.CParm table (i.e. it
% identifies which file to read/preprocess).
% parms will contain the struct stored in ns.CParm. How this parms
% struct is used is completely up to the read function.(see the example read functions below).
% Before calling the read function, ns.C checks that the file exists.
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
%           info field of the CChannel table. For example, this could contain the
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
%  To use these, add a row to CParm with a line like thise
% insert(ns.CParm,struct('ctag','eeg','fun','ephys.intan.read','description'','EEG
%                 data', 'parms',parmsStruct))
% where parmsStruct is a struct that defines what kind of preprocessing the
% read function should do before adding the data to the ns.C table.
%
% The align() function of the ns.C class is used in the analysis to
% extract trial-based data, or to align to specific events in the
% experiment.
%
% BK - June 2023

classdef C< dj.Computed
    properties (Dependent)
        keySource
    end

    methods

        function v = get.keySource(~)
            % Restricted to files with the extenstion specified in CParm
            % and the include/exclude specs in CParm.
            % This seems cumbersome, but I coudl not get a simpler join to work
            allTpl = [];
            for thisPrm= fetch(ns.CParm,'extension','include','exclude')'
                % Loop over the rows in CParm
                restrict  =struct('extension',thisPrm.extension);
                tbl = ns.File & restrict;
                if ~isempty(thisPrm.include)
                    inc = strsplit(thisPrm.include,',');
                    for i=1:numel(inc)
                        tbl = tbl & ['filename LIKE ''' inc{i} ''''];
                    end
                end

                if ~isempty(thisPrm.exclude)
                    exc = strsplit(thisPrm.exclude,',');
                    for i=1:numel(exc)
                        tbl = tbl & ['filename NOT LIKE ''' exc{i} ''''];
                    end
                end
                % Table for one row in CParm
                tbl = tbl*proj(ns.CParm&ns.stripToPrimary(ns.CParm,thisPrm));
                % Would like to concatenate this tbl with the next row but
                % this does not work with the | operator. Instead, concatenate
                % tuples of primary keys
                thisTpl = fetch(tbl);
                if isempty(allTpl)
                    allTpl = thisTpl;
                else
                    allTpl  = catstruct(1,allTpl,thisTpl);
                end
            end
            % And then restrict the full table by the set of found tuples.
            v = (ns.File*proj(ns.CParm,'fun','description','parms')) & allTpl;
        end
    end
    methods (Access=public)
        function [n,T] = toPreprocess(tbl,tag)
            % Return the number of files that stll need to be processed. If
            % a second output is requested, also returns a table with the
            % items that are on the list to be processed.
            total = tbl.keySource & struct('ctag',tag);
            done = tbl & struct('ctag',tag);
            n = count(total-done);
            if nargout>1
                T =(total-done);
            end
        end
        function [varargout] = align(tbl,pv)
            % Function to retrieve trial-based and aligned preprocessed signals.
            %
            % A sample is assigned to a trial if it occurs after the first monitor frame
            % in the trial and before the first monitor frame of the next trial.
            % (In other words the ITI is included at the *end* of each trial).
            % INPUT
            % tbl  - ns.C table
            %
            % Optional Parameter/Value pairs
            % channel  - The subset of channels ([]; means all). A cellstr
            % of channel names is valid too (or a single channel name as
            % string or char)
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
            % fetchOptions - Options passed to the fectch(ns.C)
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
            % T= align(ns.C & 'ctag=''eeg''',channel=2,align=
            % stimulusOnsetTime, start =-250,stop=1000,step=1)
            % will return a timetable with the EEG samples of channel 2,
            % from 250 ms before until 1000 ms after the stimulus onset
            % time in each trial. (stimulusOnsetTime has to be a vector
            % with a length equal to the number of trials in which each
            % element represent the onset time in ms such as returned by ns.Experiment.get()
            arguments
                tbl  (1,1) ns.C {mustHaveRows(tbl,1)}
                pv.fetchOptions {mustBeText} = ''
                pv.channel (1,:)  =[]   % The channel number (as a vector of double) or name (cellstr or string, char)
                pv.trial (1,:) double = []
                pv.start (1,1) double = 0
                pv.stop  (1,1) double = 1000
                pv.step (1,1) double  = 1;
                pv.interpolation {mustBeText} = 'linear'
                pv.crossTrial (1,1) logical = false;
                pv.align (1,:) double = []
                pv.removeArtifacts (1,1) = true
            end
           
            %% Collect data from dbase and check consistency
            sampleTime = fetch1(tbl ,'time');
            if numel(sampleTime)==3
                sampleTime = linspace(sampleTime(1),sampleTime(2),sampleTime(3))';
            end
            nrAllTrials = fetch1(ns.Experiment*tbl,'nrtrials');
            allTrials = 1:nrAllTrials;
            if isempty(pv.trial)
                trials = allTrials; % All trials
            else
                trials = intersect(pv.trial,allTrials);
            end
            
            
            % trial time zero in NS for each trial
            trialStartTime = get(ns.Experiment & tbl, 'cic','prm','firstFrame','atTrialTime',inf,'what','clocktime');

            if ~pv.crossTrial && isinf(pv.stop)
                % use maximum trial duration to setup the T table            
                pv.stop = max(diff(trialStartTime(trials)));
            end
            if isempty(pv.align)
                pv.align = zeros(size(trials)); % Align to first frame
            elseif numel(pv.align)==1 % Singleton expansion
                pv.align = repmat(pv.align,[1 numel(trials)]);
            else
                assert(numel(pv.align)==numel(trials),'Each trial should have an align time ')
            end
            % Retrieve the signal for the requested channels in the entire experiment
            % [nrSamples nrTrials]
            if isempty(pv.channel)
                channelRestriction = struct([]);
            elseif isnumeric(pv.channel)
                channelRestriction = struct('channel',num2cell(pv.channel(:))');
            elseif ischar(pv.channel) || isstring(pv.channel) ||iscellstr(pv.channel)
                channelRestriction = struct('name',pv.channel);
            elseif isstruct(pv.channel)
                channelRestriction = pv.channel;
            end
            tblChannel =ns.CChannel & proj(tbl) & channelRestriction;
            if ~isempty(pv.fetchOptions)
                channelTpl = fetch(tblChannel,'signal',pv.fetchOptions);
            else
                channelTpl = fetch(tblChannel,'signal');
            end
            signal =double([channelTpl.signal]); % Signal as matrix
            [nrSamples,nrChannels] = size(signal);

            %% Artifact removal
            if pv.removeArtifacts
                % Correction that applies to all channels
                 aTbl = ns.Artifact&tbl;
                 if exists(aTbl)
                    exptArtifacts= fetch(aTbl,'trial','start','stop');          
                    for tr=exptArtifacts.trial
                        from = sampleTime>=trialStartTime(tr);
                        if tr<nrAllTrials
                            to = sampleTime<=trialStartTime(tr+1);
                        else
                            to = from;
                        end
                        signal(from & to,:)=NaN;
                    end
                    isArtifact = any(arrayfun(@(a,b) (sampleTime>=a & sampleTime<=b),exptArtifacts.start,exptArtifacts.stop,'UniformOutput',true),2);
                    signal(isArtifact,:) = NaN;
                 end
                
                % Artifacts found in individual channels
                 acTbl =  (ns.ArtifactChannel & proj(aTbl ))& struct('channel',{channelTpl.channel});
                 if exists(acTbl)
                    channelArtifacts= fetch(acTbl,'trial','start','stop');     
                    for ch= 1:nrChannels
                        for tr=channelArtifacts(ch).trial
                            from = sampleTime>=trialStartTime(tr);
                            if tr<nrAllTrials
                                to = sampleTime<=trialStartTime(tr+1);
                            else
                                to = from;
                            end
                        signal(from & to,ch)=NaN;
                        end
                        isArtifact = any(arrayfun(@(a,b) (sampleTime>=a & sampleTime<=b),channelArtifacts(ch).start,channelArtifacts(ch).stop,'UniformOutput',true),2);
                        signal(isArtifact,ch) = NaN;
                    end                                        
                 end            
            end
            % Fetch the preprocessing tag. The dimensions in the timetable
            % will be named after this.
            tag = fetch1(ns.CParm & tbl,'ctag');
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
                T.Properties.DimensionNames = {'TrialTime',tag};
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
                error('File %s does not exist. Cannot create C',filename);
            end

            % Get the specified preprocessing parms
            prepParms = fetch(ns.CParm & key,'*');
            % Call the prep function
            [signal,time,channelInfo,recordingInfo] = feval(prepParms.fun,key,prepParms.parms);
            [nrSamples,nrChannels] = size(signal);

            assert(nrSamples==numel(time) || (numel(time)==3 &&time(3)==nrSamples),'The number of rows in the preprocessed signal does not match the number of time points ')
            assert(nrChannels==numel(channelInfo),'The number of columns in the preprocessed signal does not match the number of channels')

            channels =[channelInfo.nr];

            % Create tuples and insert.
            tpl = mergestruct(key,...
                struct('time',time, ...
                'nrsamples',nrSamples,...
                'info',recordingInfo));
            insert(tbl,tpl)

            % Create tpls for each of the channels and insert
            channelsTpl = mergestruct(key,...
                struct('signal',num2cell(single(signal),1)',...
                'channel',num2cell(channels(:))));

            if ~isempty(channelInfo)
                for i=1:numel(channelInfo)
                    channelsTpl(i).channelinfo = channelInfo(i);
                    if isfield(channelInfo,'name')
                        channelsTpl(i).name = channelInfo(i).name;
                    end
                end
            end

            % Chunking the inserts to avoid overloading the server
            chunkSize = 128; % This should probably be user configurable (e.g., NS_MAXUPLOAD)
            tic;
            fprintf('Uploading to server ')
            for i=1:chunkSize:nrChannels
                fprintf('.')
                thisChunk = i:min(nrChannels,i+chunkSize-1);
                insert(ns.CChannel,channelsTpl(thisChunk));
            end
            fprintf('Done in %d seconds.\n.',round(toc))

        end
    end

    methods (Access=public)

        function out = plot(cTbl,channel, expt, grouping,pv)
            % Plot time courses, spectra.
            %
            % Each roi in the table will be shown as a
            % separate tile, each condition a line in the plot.  Time
            % courses are scaled to the 99th percentile across all
            % responses and de-meaned per condition. Hence, the mean
            % response is lost, but the relative response modulation in
            % each condition is maintained.
            %
            % tbl - ns.C*ns.CChannel table
            % expt - A single ns.Experiment (tuple or table)
            % condition - Specify how trials should be grouped into conditions:
            %               []  - Pool over all trials
            %               A ns.Condition table - pool per condition
            %               A vector of trials - Pool over only these trials
            %               A cell array with vectors of trials. Pool over
            %               each set of trials
            %               ns.Dimension  - group by the conditions in the
            %               dimension
            %
            % 'fun' - By default this function visualizes the mean across trials
            %       together with shading reflecting the standard error. To use something else,
            %       pass a function that, when passed a matrix with [nrTimePoints nrTrials] , returns one
            %       value and an error bar for each row.
            %
            % 'name'  Name of the conditions
            % 'start' - Start time in seconds
            % 'step'  - Step time in seconds
            % 'stop' - Stop time in seconds
            % 'interpolation' - Interpolation method ['linear']
            % 'crossTrial ' - Allow start/stop to cross to the
            %               previous/next trial.
            % 'mode'  ["TIMECOURSE"], TUNING, SPECTRUM ,"RASTER"
            % 'perTrial'  -Show individual trials [false]
            % 'prctileMax'  Percentile that is used to scale responses for
            %               visualization [95]
            % Spectrum Options
            % 'evoked' Set  to true to show evoked power instead of total
            % power.

            arguments
                cTbl (1,1) ns.C {mustHaveRows}
                channel   % A ns.CChannel or a CChannel based restriction
                expt (1,1)
                grouping = []
                pv.removeArtifacts (1,1) = true
                pv.trial = [] 
                pv.fun (1,1) = @(x)(deal(mean(x,2,"omitnan"),std(x,0,2,"omitnan")./sqrt(sum(~isnan(x),2))));
                pv.name {mustBeText} = num2str(1:numel(grouping))
                pv.start (1,1) double = 0
                pv.stop (1,:) double =  2000
                pv.step  (1,1) double = 100;
                pv.align (1,:) double = []
                pv.interpolation {mustBeText} = 'linear';
                pv.modality {mustBeText} = 'spikes';
                pv.averageOverChannels (1,1)  logical = false;
                pv.mode (1,:) {mustBeMember(pv.mode,["COHERENCE", "RASTER", "TIMECOURSE","EVOKED","TOTAL"])} = "TIMECOURSE"
                pv.crossTrial (1,1) logical = true;
                pv.fetchOptions {mustBeText} = ''
                pv.perTrial (1,1) logical = false;
                pv.prctileMax (1,1) double {mustBeInRange(pv.prctileMax,0,100)} = 95;
                % Layout
                pv.compact = false;
                pv.fig  = []; % Pass a handle to a figure to reuse/add
                % Spectrum options
                pv.evoked (1,1) logical = false;
                pv.options cell = {}; % Cell array of parameter value pairs passed to pspectrum
            end
            exptTpl = fetch(ns.Experiment & expt);
            % Convert the channel specification to a restriction
            if isempty(channel)
                channelRestriction = struct([]);
            elseif isnumeric(channel)
                channelRestriction = struct('channel',num2cell(channel(:))');
            elseif ischar(channel) || isstring(channel) ||iscellstr(channel)
                channelRestriction = struct('name',channel);
            elseif isstruct(channel)
                channelRestriction = channel;
            end

            % Determine how to group the trials
            if iscell(grouping)
                % A cell array of trial numbers with user-specified
                % groupins
                trials = grouping;
                names = pv.name;
                conditionOrder = 1:numel(grouping);
            else
                % grouping identifies a dimension
                conditions = ns.Dimension*ns.DimensionCondition & grouping & exptTpl;
                [trials,names,conditionX] = fetchn(conditions,'trials','name','value');
                conditionX = [conditionX{:}];
                conditionX = cell2mat(conditionX);
                [~,conditionOrder] =sort(conditionX);
            end

            if ~isempty(pv.trial)
                trials = cellfun(@(x) intersect(x,pv.trial),trials,'uni',false);
            end
            nrConditions = numel(trials);
            nrChannels = count(cTbl*(ns.CChannel & channelRestriction) & expt);

            % Loop over Channels
            perTrial =cell(nrConditions,nrChannels);  % Collect per trial to return
            for channelCntr = 1:nrChannels
                % First, read the data for each condition and store this
                % locally. Because this takes time and the same data may be
                % used for multiple "mode"s, this saves time.
                for c= conditionOrder
                    % Loop over conditions in the order specified above
                    if pv.averageOverChannels || any(contains(pv.mode,"COHERENCE"))
                        % Fetch all channels at once
                        fetchOptions = pv.fetchOptions;
                    else
                        % Pick the current channel (i.e a row in
                        % cTbl&expt&channelRestriction  , use pv.fetchOptions
                        % to force a sort order)
                        fetchOptions = ['LIMIT 1 OFFSET '  num2str(channelCntr-1) ' ' pv.fetchOptions];
                    end
                     if isempty(pv.align)
                            alignTime = zeros(size(trials{c})); % Align to first frame
                     elseif numel(pv.align)==1 % Singleton expansion
                            alignTime = repmat(pv.align,[1 numel(trials{c})]);
                     else
                            alignTime = pv.align(trials{c});
                    end


                    [y,time] = align(cTbl&exptTpl,removeArtifacts = pv.removeArtifacts, channel =channelRestriction, trial=trials{c},fetchOptions= fetchOptions,align = alignTime, crossTrial =pv.crossTrial,start=pv.start,stop=pv.stop,step=pv.step,interpolation =pv.interpolation);
                    if pv.averageOverChannels
                        % Average
                        y = mean(y,2,"omitnan"); % Average over rois
                    end
                    if isempty(y);continue;end
                    % Store the per-trial data (used only by some modes)
                    perTrial{c,channelCntr} = y;
                end

                % Now analyze the data for this channel
                for mode= pv.mode
                    % Initialize the y-values, errors and time variables
                    m = [];  % Average over trial
                    e = [];
                    allX = [];
                    if isempty(pv.fig)
                        figByName(sprintf('%s (Ch#%d) -%s on %s@%s',mode,channelCntr,exptTpl.subject ,exptTpl.session_date ,exptTpl.starttime));
                        clf;
                        layout = tiledlayout('flow');
                        if pv.compact
                            layout.Padding ="tight";
                        end
                    else
                        % Assume that this function has been called before and the
                        % layout has already been created.
                        if isempty(pv.fig.Children)
                            layout = tiledlayout('flow');
                        else
                            layout = pv.fig.Children;
                        end
                        figure(pv.fig);
                    end
                    
                    for c= conditionOrder
                        if c==1
                            if isempty(pv.fig) || numel(layout.Children)<channelCntr
                                nexttile;
                            else
                                axes(layout.Children(channelCntr));
                            end
                        end
                        y = perTrial{c,channelCntr};

                        % Post-process depending on mode
                        switch upper(mode)
                            case {"TOTAL","EVOKED"}
                                % TOTAL or EVOKED Power
                                y = y-mean(y,1,"omitnan");
                                if upper(pv.mode) =="EVOKED"
                                    y = mean(y,2,"omitnan");
                                end
                                y(isnan(y)) =0;
                                [pwr,freq] = pspectrum(y,time,'power',pv.options{:});
                                % [ft,freq] = fftReal(y,1./seconds(time(2)-time(1)));
                                % pwr = ft.*conj(ft);
                                [thisM,thisE] = pv.fun(pwr); % Average over trials
                                thisX = freq;
                                thisM = thisM.*freq;
                            case {"TIMECOURSE", "RASTER"}
                                % Average over trials  in the condition
                                [thisM,thisE] = pv.fun(y);
                                thisX = time;
                            case "COHERENCE"
                                % TODO Determine coherence across channels
                                y = y- mean(y,1,"omitnan"); % Remove mean
                                y(isnan(y)) = 0; % Remove nans
                                for tr =  1:size(y,2)
                                    [thisC(:,:,:,tr),phi,S12,freq] = cohmatrixc(squeeze(y(:,tr,:)),struct('tapers',[3 5],'pad',0,'Fs',1./pv.step));
                                end
                                thisM = mean(thisC,4);
                                thisX = time; % not sure yet
                        end
                        m = catpad(m,thisM);  % Cat as next column allow different rows (padded with NaN at the end)
                        e  =catpad(e,thisE);
                        if numel(thisX)>numel(allX)
                            allX =thisX;
                        end
                    end


                    %% Visualize per channel
                    nrX = numel(allX);
                    switch upper(mode)
                        case {"TOTAL","EVOKED"}
                            grandMax = prctile(abs(m(:)),pv.prctileMax );
                            grandMin = prctile(abs(m(:)),100-pv.prctileMax );
                            m = (m-grandMin)./(grandMax-grandMin);
                            e = e./(grandMax-grandMin);
                            % Add the conditionNr so that each m column has a mean of
                            % conditionNr and can be plotted on the same axis, with
                            % conditions discplaced vertically from each other.
                            m = m + repmat(1:nrConditions,[nrX 1]);
                            [h,hErr] = ploterr(allX,m,e,'linewidth',2,'ShadingAlpha',0.5);
                            hold on
                            % Show "zero" line
                            hh = plot(allX,repmat(1:nrConditions,[nrX 1]),'LineWidth',0.5);
                            [hh.Color] =deal(h.Color);
                            ylim([0 nrConditions+1])
                            set(gca,'yTick',1:nrConditions,'yTickLabel',names(conditionOrder))
                            xlabel 'Frequency (Hz)'
                            h =legend(h,names(conditionOrder));
                            h.Interpreter  = 'None';
                        case "RASTER"
                            % Using an image to show the time course.
                            grandMax = max(cellfun(@(x) prctile(x(:),pv.prctileMax ),perTrial(:,channelCntr)));
                            grandMin = min(cellfun(@(x) prctile(x(:),100-pv.prctileMax ),perTrial(:,channelCntr)));
                            nrX = numel(allX);
                            cmap = [hot(255);0 0 1];
                            I =[];
                            for c=1:nrConditions
                                I = cat(2,I,perTrial{c,channelCntr},nan(nrX,1));
                            end
                            I = ((I-grandMin)./(grandMax-grandMin))';
                            % Clamp
                            I(I<0) = 0;
                            I(I>1) = 1;
                            I= round(I*255);
                            I(isnan(I))=256;
                            nrTrials = size(I,1);
                            image(seconds(allX),1:nrTrials, I,'CDataMapping','direct');
                            colormap(cmap)
                            nrTrialsPerCondition = cellfun(@(x) size(x,2),perTrial(:,channelCntr))+1;
                            leftEdge = [0; cumsum(nrTrialsPerCondition(1:end-1))];
                            middleOfCondition = leftEdge+nrTrialsPerCondition./2;
                            set(gca,'yTick',middleOfCondition,'yTickLabel',names(conditionOrder))
                            if pv.compact
                                set(gca,'XTick',[]);
                            else
                                xlabel 'Time (s)'
                                ylabel('Conditions')
                                h = colorbar;
                                set(h,'YTick',0:50:250,'YTickLabel',round(grandMin +(0:50:250)*(grandMax-grandMin)/255) )
                                ylabel(h,'Response')
                            end
                        case "TIMECOURSE"
                            %% One time series line perCondition , spaced vertically.
                            % Scale each condition to the grandMax
                            if pv.perTrial
                                grandMax = max(cellfun(@(x) prctile(x(:),pv.prctileMax ),perTrial(:,channelCntr)));
                                grandMin = min(cellfun(@(x) prctile(x(:),100-pv.prctileMax ),perTrial(:,channelCntr)));
                            else
                                grandMax = prctile(abs(m(:)),pv.prctileMax );
                                grandMin = prctile(abs(m(:)),100-pv.prctileMax );
                            end
                            m = (m-grandMin)./(grandMax-grandMin);
                            e = e./(grandMax-grandMin);
                            % Add the conditionNr so that each m column has a mean of
                            % conditionNr and can be plotted on the same axis, with
                            % conditions discplaced vertically from each other.
                            m = m + repmat(1:nrConditions,[nrX 1]);
                            [h,hErr] = ploterr(allX,m,e,'linewidth',2,'ShadingAlpha',0.5);
                            hold on
                            % Show "zero" line
                            hh = plot(allX,repmat(1:nrConditions,[nrX 1]),'LineWidth',0.5);
                            [hh.Color] =deal(h.Color);
                            ylim([0 nrConditions+1])
                            set(gca,'yTick',1:nrConditions,'yTickLabel',names(conditionOrder))
                            xlabel 'Time (s)'
                            ylabel 'Response per condition'
                            if pv.perTrial
                                [hh.Color] = deal([0 0 0]);
                                [h.Color] = deal([0 0 0]);
                                [hErr.FaceColor] = deal([0 0 0]);
                                colorOrder = get(gca,'ColorOrder');
                                for c=1:nrConditions
                                    plot(allX,c+ perTrial{:,channelCntr}./(grandMax-grandMin),'Color',colorOrder(mod(c-1,size(colorOrder,1))+1,:),'LineWidth',.5)
                                end
                            end
                    end
                end
                if pv.averageOverChannels || mode=="COHERENCE"
                    break;
                end
            end
            if nargout>0
                out =perTrial;
            end
        end

        function [hdr,data,evts] = fieldtrip(tbl,pv)
            arguments
                tbl (1,1) ns.C
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
                tpl = fetch(tbl*(proj(ns.CChannel,'info->channelInfo','signal') &  key) & struct('channel',num2cell(pv.channel)'),'*') ;
                nrChannels = numel(tpl);
                channelInfo = [tpl.channelInfo];
                nrTrials = fetch1(ns.Experiment &key,'nrtrials');
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