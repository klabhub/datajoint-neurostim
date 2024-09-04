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
        function [time,trial,value] = eventTrialTime(tbl,pv)
            arguments
                tbl (1,1) ns.C
                pv.mode (1,1) string {mustBeMember(pv.mode,["AFTER" "BEFORE" "NEAREST"])} = "AFTER"
                pv.value  (1,:) double = NaN
            end
            for tpl = fetch(tbl*ns.CChannel,'time','nrsamples','signal')'
                zeroTime = get(ns.Experiment & tpl,'cic','prm','firstFrame','what','clocktime');
                if ~isnan(pv.value(1))
                    stay = ismember(tpl.signal,pv.value);
                else
                    stay = true(size(tpl.time));
                end
                dt = tpl.time(stay)- zeroTime';
                modDt = dt;
                if pv.mode =="NEAREST"
                    modDt = abs(modDt);
                elseif pv.mode == "AFTER"
                    modDt(modDt<0)=inf;
                elseif pv.mode =="BEFORE"
                    modDt(modDt>0)=inf;
                end
                % For each event find the closest firstframe eventtime
                [~,trial] = min(modDt,[],2);
                ix = sub2ind(size(dt),(1:sum(stay))',trial);
                time = dt(ix);
                value =tpl.signal(stay);
            end

        end

        function v = get.keySource(~)
            % Restricted to files with the extenstion specified in CParm
            % and the include/exclude specs in CParm.
            % This seems cumbersome, but I coudl not get a simpler join to work
            allTpl = [];
            for thisPrm= fetch(ns.CParm,'extension','include','exclude')'
                % Loop over the rows in CParm
                restrict  =struct('extension',thisPrm.extension);                
                tbl = ns.File & restrict & analyze(ns.Experiment,strict=false); % Only files in experiments that should be analyzed
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
        function [t,dt] = sampleTime(tbl)
            % Determine time and time step of this ns.C entry.
            arguments
                tbl (1,1) {mustHaveRows(tbl,1)}
            end
            t= double(fetch1(tbl ,'time'));
            if numel(t)==3
                t= linspace(t(1),t(2),t(3))';
            end
            dt = mode(diff(t));
        end

        function varargout = plot(cTbl,pv)
            % function out = plot(cTbl,pv)
            % Plot time courses, spectra.
            %
            % Each channel will be shown as a
            % separate tile, each condition a line in the plot.  Time
            % courses are scaled to the 99th percentile across all
            % responses and de-meaned per condition. Hence, the mean
            % response is lost, but the relative response modulation in
            % each condition is maintained.
            %
            % cTbl - ns.C  table
            % channel - Which channel to show. This can be a number to
            % match with CChannel.channel, a string to match with
            % CChannel.name or any restriction on the CChannel table.
            %
            %
            % grouping- Specify how trials should be grouped into conditions:
            %               []  - All trials are considered a single
            %               condition (Default).
            %               A cell array of trial numbers identifies trials
            %               as belonging to specific conditions, each will
            %               be plotted separately.
            %               A string idenfities a dimension  in
            %               ns.Dimension; the different conditions in this
            %               dimension will be shown separately.
            %               A struct identifies a restriction on
            %               ns.Dimension.
            %
            % 'removeArtifacts' [true] - Use the ns.Artifacts table to
            % select trials and time periods without artifacts, and plot
            % only those.
            %
            % 'average' - The function used to determine an average across trials. By default this
            %           determines the mean and standard errors.
            %       together with shading reflecting the standard error. To use something else,
            %       pass a function that, when passed a matrix with
            %       [nrTimePoints nrTrials] , returns two column vector outputs; the
            %       average value and an error.
            %
            % 'name'  Name of the conditions to be used when grouping
            %           does not use the ns.Dimension table.
            % 'start' - Start time in milliseconds
            % 'step'  - Step time in milliseconds
            % 'stop' - Stop time in milliseconds
            % 'interpolation' - Interpolation method ['nearest']. See
            %                   ns.C/align
            % 'crossTrial ' - Allow start/stop to cross to the
            %               previous/next trial.
            % 'averageOverChannels [false] - Average over channels, or show
            %                   each channel in a separate tile.
            % 'fetchOptions' - Fetch options to use when fetching the data
            %                   from ns.C
            % 'mode'  ["TIMECOURSE"], RASTER, TOTAL, EVOKED, COHERENCE
            % 'prctileMax'  Percentile that is used to scale responses for
            %               visualization [95]
            % 'padding'  - Padding to use in the tiledlayout  ['tight']
            % 'forceFig'  - Force creating a new figure [false]
            % Spectrum Options
            % 'pspectrum' A cell array of options passed to pspectrum in
            % TOTAL or EVOKED modes
            %
            % EXAMPLE:
            % Determine the onset of a visual stimulus in experiment expt
            % onset  = get(ns.Experiment & expt,'flicker','prm','startTime','what','trialtime','atTrialTime',inf);
            % Plot the timecourse of channels 1 and 2 of the eeg C data, grouped by the
            % isFlick dimension, aligned to the stimulus onset:
            % plot(ns.C & 'ctag=''eeg''',[1 2 ],expt,"isFlick",prctileMax=99, start =-100,stop =500,align=onset,mode ="TIMECOURSE",perTrial =true);

            arguments
                cTbl (1,1) ns.C {mustHaveRows}
                pv.channel   = []  % A ns.CChannel or a CChannel based restriction                
                pv.grouping = []
                pv.groupingName = {};
                pv.removeArtifacts (1,1) = true
                pv.trial = []
                pv.fun = @(x)(x) % A function to apply to the data obtained from ns.C
                pv.average (1,1) = @(x)(deal(mean(x,2,"omitnan"),std(x,0,2,"omitnan")./sqrt(sum(~isnan(x),2))));
                pv.name {mustBeText} = ""
                pv.start (1,1) double = 0
                pv.stop (1,:) double =  inf
                pv.step  (1,1) double = 0;
                pv.align (1,:) double = []
                pv.interpolation {mustBeText} = 'nearest';

                pv.averageOverChannels (1,1)  logical = false;
                pv.mode (1,:) {mustBeMember(pv.mode,["COHERENCE", "RASTER", "TIMECOURSE","EVOKED","TOTAL"])} = "TIMECOURSE"
                pv.crossTrial (1,1) logical = false;
                pv.fetchOptions {mustBeText} = ''
                pv.prctileMax (1,1) double {mustBeInRange(pv.prctileMax,0,100)} = 95;
                % Layout
                pv.padding string = "compact";
                pv.forceFig  = false
                pv.useImage = false;
                % Spectrum options
                pv.pspectrum cell = {}; % Cell array of parameter value pairs passed to pspectrum
            end
            

            [T,conditionName,channelNr] = align(cTbl,channel=pv.channel,grouping=pv.grouping,trial=pv.trial, ...
                            start=pv.start,stop=pv.stop,step=pv.step, interpolation = pv.interpolation,crossTrial =pv.crossTrial,...
                            average=pv.average,averageOverChannels= pv.averageOverChannels, ...
                            fun=pv.fun,removeArtifacts =pv.removeArtifacts ,...
                            align = pv.align,fetchOptions = pv.fetchOptions);
            nrChannels = numel(channelNr);
            if isa(T,"timetable");T={T};end
            out= cellfun(@isempty,T);
            T(out) = [];
            conditionName(out)=[];
            if ~isempty(pv.groupingName)
                conditionName =pv.groupingName;
            end
                

            nrConditions= size(conditionName,1);
            exptTpl =fetch(ns.Experiment &cTbl);
            ctag = fetch1(cTbl,'ctag');
            for channelCntr = 1:nrChannels
                thisChannelName= string(channelNr(channelCntr));
                
                % Now analyze the data for this channel
                for mode= pv.mode
                    % Initialize the y-values, errors and time variables
                    m = [];  % Average over trial
                    e = [];
                    allX = [];
                    if pv.forceFig
                        uid = "uid: " + string(randi(1e10));
                    else
                        uid = "";
                    end
                    hFig = figByName(sprintf('%s (%s) -S:%s on %s@%s %s',mode,ctag,exptTpl.subject ,exptTpl.session_date ,exptTpl.starttime,uid));
                    if channelCntr==1
                        clf;
                        layout = tiledlayout('flow');
                        layout.Padding =pv.padding;
                    else
                        layout  = hFig.Children;
                    end
                    hAx = findobj(layout,'Type','axes');                    
                    if channelCntr > numel(hAx)
                        nexttile;
                    else
                        axes(hAx(1)) %#ok<LAXES>
                    end
                    nrTrialsPerCondition =nan(1,nrConditions);
                    for c= 1:nrConditions                        
                        [y,time] = timetableToDouble(T{c});
                        % Post-process depending on mode
                        switch upper(mode)
                            case {"TOTAL","EVOKED"}
                                % TOTAL or EVOKED Power
                                y = y-mean(y,1,"omitnan");
                                if upper(pv.mode) =="EVOKED"
                                    y = mean(y,2,"omitnan");
                                end
                                y(isnan(y)) =0;
                                % fb = cwtfilterbank("SamplingFrequency",1000./pv.step,"SignalLength",size(y,1))
                                % for tr = 1:size(y,2)
                                %     s(:,:,tr) =cwt(y(:,tr),fb);
                                % end

                                [pwr,freq] = pspectrum(y,time,'power',pv.pspectrum{:});
                                % [ft,freq] = fftReal(y,1./seconds(time(2)-time(1)));
                                % pwr = ft.*conj(ft);
                                [thisM,thisE] = pv.average(pwr); % Average over trials
                                thisX = freq;
                                thisM = thisM.*freq;
                            case "TIMECOURSE"
                                % Average over trials  in the condition
                                [thisM,thisE] = pv.average(y);
                                thisX = time;
                            case "RASTER"
                                thisM = y;
                                nrTrialsPerCondition(c) = size(y,2);
                                thisE =  [];
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
                            if pv.useImage
                                imagesc(allX,1:numel(x),m')
                                set(gca,'yTick',1:numel(x),'yTickLabel',num2str(x'))
                                xlabel 'Frequency (Hz)'
                                axis xy
                                hC = colorbar;
                                ylabel(hC,sprintf('Power (%.2f-%.2f)',grandMin,grandMax))
                                hC.YTick = [];
                            else
                                % Use line plots per condition
                                % Add the conditionNr so that each m column has a mean of
                                % conditionNr and can be plotted on the same axis, with
                                % conditions discplaced vertically from each other.
                                m = m + repmat(1:nrConditions,[nrX 1]);
                                [h,~] = ploterr(allX,m,e,'linewidth',2,'ShadingAlpha',0.5);
                                hold on
                                % Show "zero" line
                                hh = plot(allX,repmat(1:nrConditions,[nrX 1]),'LineWidth',0.5);
                                [hh.Color] =deal(h.Color);
                                ylim([0 nrConditions+1])
                                set(gca,'yTick',1:nrConditions,'yTickLabel',conditionName)
                                xlabel 'Frequency (Hz)'                                
                            end
                        case "RASTER"
                            % Using an image to show the time course.
                            grandMax = prctile(abs(m(:)),pv.prctileMax );
                            grandMin = prctile(abs(m(:)),100-pv.prctileMax );                          
                            nrX = numel(allX);
                            cmap = [hot(255);0 0 1];
                            I =[];
                            trialsSoFar = 0;
                            for c=1:nrConditions
                                I = cat(2,I,m(:,trialsSoFar+(1:nrTrialsPerCondition(c))),nan(nrX,1));
                                trialsSoFar = trialsSoFar+ nrTrialsPerCondition(c);
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
                            axis xy
                            nrLinesPerCondition = nrTrialsPerCondition+1; % Marker to separate conditions
                            leftEdge = [0 cumsum(nrLinesPerCondition(1:end-1))];
                            middleOfCondition = leftEdge+nrLinesPerCondition./2;
                            set(gca,'yTick',middleOfCondition,'yTickLabel',conditionName)
                            if pv.padding=="tight"
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
                            set(gca,'yTick',1:nrConditions,'yTickLabel',conditionName)
                            xlabel 'Time (s)'
                            ylabel 'Response per condition'

                    end
                end
                if pv.averageOverChannels || mode=="COHERENCE"
                    break;
                end
                title(sprintf('Channel: %s',thisChannelName))
            end
            if nargout>0
                varargout{1} =m;
                if nargout>1
                    varargout{2} = allX;                    
                    if nargout >2
                        varargout{3} = time;     
                        if nargout >3
                            varargout{4} =conditionName;
                        end
                    end
                end
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

        function [T,conditionValue,channelNr] = align(tbl,pv)
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
            % start  - Extract signals startint at this time (ms) [0]
            % stop   - Last time point to extract. (ms) [inf; end of trial]
            % step   - Step size in seconds. (ms)   [To use the native
            % resolution of the ns.C row, use 0, to average or interpolate,
            % specify a time].
            % align - The time in each trial that is considered 0. By
            % default this is the time at which the first monitor frame
            % became visible. By specifying the time of an event (e.g.,
            % stimulus onset), the data can be aligned to that event.
            %
            % interpolation -  enum('nearest','linear','spline','pchip','makima')
            %               Interpolation method; see timetable/synchronize. ['nearest']
            %
            % crossTrial - Allow values to be returned that are from the
            %               trials before or after. When this is false,
            %               only times/samples between the firstFrame event in the trial and
            %               the firstFrame event of the next trial will be
            %               returned.
            % fetchOptions - Options passed to the fectch(ns.C)
            %               call.(For instance 'LIMIT 1')
            %
            % OUTPUT
            % T     = timetable with each column a trial. Time is in seconds
            %           relative to the first frame of the trial.
            %          Channels are along the columns of the
            %           elements of the table. 
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
                pv.channel  =[]   %
                pv.grouping = []
                pv.averageOverChannels (1,1) logical =false
                pv.fun (1,1) function_handle = @(x)(x)  % A function to apply to the data
                pv.trial (1,:) double = []
                pv.start (1,1) double = 0
                pv.stop  (1,1) double = inf
                pv.step (1,1) double  = 0;
                pv.interpolation {mustBeText} = 'nearest'
                pv.crossTrial (1,1) logical = false;
                pv.align (1,:) double = []
                pv.removeArtifacts (1,1) = true
            end

            % Expt info.
            exptTpl = fetch(ns.Experiment & tbl,'nrtrials');
            trialStartTime = get(ns.Experiment & tbl, 'cic','prm','firstFrame','atTrialTime',inf,'what','clocktime');
            [t,dt] = sampleTime(tbl);
            if pv.step==0
                pv.step = dt;
            end


            %% Convert the channel specification to a restriction
            if isempty(pv.channel)
                channelRestriction = struct([]);
            elseif isnumeric(pv.channel)
                channelRestriction = struct('channel',num2cell(pv.channel(:))');
            elseif ischar(pv.channel) || isstring(pv.channel) ||iscellstr(pv.channel)
                channelRestriction = struct('name',cellstr(pv.channel)');
            elseif isstruct(pv.channel)
                channelRestriction = pv.channel;
            end
            nrChannels = count(tbl*(ns.CChannel & channelRestriction) & exptTpl);
            assert(nrChannels >0,'No matching channels found. Nothing to do.\n');
            [channelNr,channelName] = fetchn((tbl&exptTpl)*(ns.CChannel&channelRestriction),'channel','name');

            %% Group and select trials
            if iscell(pv.grouping)
                % A cell array of trial numbers with user-specified
                % groupins
                trials = pv.grouping;
                conditionOrder = 1:numel(pv.grouping);
                conditionValue = (1:numel(pv.grouping))';
                conditionName = string(conditionValue);
            elseif isempty(pv.grouping)
                trials= {1:exptTpl.nrtrials};% All as one group
                conditionValue =1;
                conditionName = string(conditionValue);
                conditionOrder =1;
            else
                if isstring(pv.grouping) || ischar(pv.grouping)
                    groupingRestriction = struct('dimension',pv.grouping);
                elseif isstruct(pv.grouping)
                    groupingRestriction  = pv.grouping;
                else
                    error('Grouping must be a string or a struct')
                end
                % grouping identifies a dimension
                conditions = (ns.Dimension & groupingRestriction & exptTpl) *ns.DimensionCondition;
                if ~exists(conditions)
                    fprintf('No conditions in dimension %s. Nothing to plot.\n',groupingRestriction.dimension);
                    return
                end
                [trials,conditionName,conditionX] = fetchn(conditions,'trials','name','value');
                conditionX = cat(1,conditionX{:}); % value can have multiple columns; cat along conditions
                conditionX = cell2mat(conditionX);
                [conditionValue,conditionOrder] =sortrows(conditionX);
            end
            if ~isempty(pv.trial)
                % Limit to the specified trials
                trials = cellfun(@(x) intersect(x,pv.trial),trials,'uni',false);
            end
            nrConditions = numel(trials);
            if isinf(pv.stop)
                if ~pv.crossTrial
                    % use maximum trial duration to setup the T table
                    pv.stop = max(diff(trialStartTime));
                else
                    error('''stop'' cannot be inf when crossTrial =true')
                end
            end


            %% Initialize
            tPerCondition =cell(nrConditions,1);
            

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
                        from = t>=trialStartTime(tr)-dt;
                        if tr<nrAllTrials
                            to = t<=trialStartTime(tr+1)+dt;
                        else
                            to = from;
                        end
                        signal(from & to,:)=NaN;
                    end
                    isArtifact = any(arrayfun(@(a,b) (t>=a-dt & t<=b+dt),exptArtifacts.start,exptArtifacts.stop,'UniformOutput',true),2);
                    signal(isArtifact,:) = NaN;
                end

                % Artifacts found in individual channels
                acTbl =  (ns.ArtifactChannel & proj(aTbl ))& struct('channel',{channelTpl.channel});
                if exists(acTbl)
                    channelArtifacts= fetch(acTbl,'trial','start','stop');
                    for ch= 1:nrChannels
                        for tr=channelArtifacts(ch).trial
                            from = t>=trialStartTime(tr)-dt;
                            if tr<nrAllTrials
                                to = t<=trialStartTime(tr+1)+dt;
                            else
                                to = from;
                            end
                            signal(from & to,ch)=NaN;
                        end
                        isArtifact = any(arrayfun(@(a,b) (t>=a-dt & t<=b+dt),channelArtifacts(ch).start,channelArtifacts(ch).stop,'UniformOutput',true),2);
                        signal(isArtifact,ch) = NaN;
                    end
                end
            end

            %Apply the fun to preprocess/modify the signal (e.g. abs())
            signal  = pv.fun(signal);

          
            %% Align
            if nrSamples==0||nrChannels==0
                T= timetable; % Empty
            else

                % Read the data for each condition
                for c= conditionOrder(:)'
                    % Loop over conditions in the order specified above

                    if isempty(pv.align)
                        alignTrialTime = zeros(size(trials{c})); % Align to first frame
                    elseif numel(pv.align)==1 % Singleton expansion
                        alignTrialTime = repmat(pv.align,[1 numel(trials{c})]);
                    elseif numel(pv.align) == exptTpl.nrtrials
                        alignTrialTime = pv.align(trials{c});
                    else
                        error('align can be empty , a singleton, or a vector with times for each trial in the experiment')
                    end

                    % Setup the new time axis for the results
                    newTimes = milliseconds(pv.start:pv.step:pv.stop)';
                    nrTimes  = numel(newTimes);

                    % Create a timetable with the activity per trial
                    nrTrials = numel(trials{c});
                    varNames = "Trial" + string(trials{c});
                    T =timetable('Size',[nrTimes nrTrials],'RowTimes',newTimes,'VariableTypes',repmat("doublenan",[1 nrTrials]),'VariableNames',varNames);

                    % Loop over trials to collect the relevant samples
                    trialOut =false(1,nrTrials);   
                    alignNsTime = nan(1,nrTrials);
                    for trCntr=1:nrTrials
                        if isinf(alignTrialTime(trCntr)) || isnan(alignTrialTime(trCntr))
                            fprintf('Align even did not occur in trial %d. Skipping trial.\n',trials{c}(trCntr))
                            trialOut(trCntr) = true;
                            continue;
                        end
                        thisTrial =trials{c}(trCntr);

                        % Limits if crossTrial is allowed
                        startLimit =  trialStartTime(thisTrial) + alignTrialTime(trCntr)+pv.start;
                        stopLimit =  trialStartTime(thisTrial) + alignTrialTime(trCntr)+pv.stop;
                        if ~pv.crossTrial
                            startLimit =  max(trialStartTime(thisTrial),startLimit); % No samples before trialStartTime(thisTrial)
                            if thisTrial<exptTpl.nrtrials
                                stopLimit = min(trialStartTime(thisTrial+1),stopLimit); % No sample after the start of next trial (last trial includes everything until the end of recording).
                            end
                        end
                        staySamples = t >= startLimit & t < stopLimit;

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
                        alignNsTime(trCntr) = trialStartTime(thisTrial)+ alignTrialTime(trCntr);
                        trialTime = t(staySamples)-alignNsTime(trCntr); %
                        thisT = timetable(milliseconds(trialTime),signal(staySamples,:)); % The table for this trial, at the original sampling rate.
                        % Now retime the table to the new time axis. Never extrapolation
                        thisT = retime(thisT,newTimes,pv.interpolation,'EndValues',NaN);
                        T.(varNames(trCntr)) = table2array(thisT);
                    end
                    T(:,trialOut) = [];                   
                    T= addprop(T,["alignTime" ],repmat("variable",[1 1]));
                    T.Properties.CustomProperties.alignTime = alignNsTime;
                    tPerCondition{c} =T;
                end

            end

            if nrConditions==1
                % Only one condition requested, return the timetable.
                T = tPerCondition{1};
            else
                % Return the cell array of timetables
                T = tPerCondition;
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
            chunkSize = 1; % This should probably be user configurable (e.g., NS_MAXUPLOAD)
            tic;
            fprintf('Uploading to server ')
            for i=1:chunkSize:nrChannels
                fprintf('.')
                if mod(i,80)==0;fprintf('\n');end
                thisChunk = i:min(nrChannels,i+chunkSize-1);
                insert(ns.CChannel,channelsTpl(thisChunk));
            end
            fprintf('Done in %d seconds.\n.',round(toc))

        end
    end




end