classdef (Abstract) cache < handle
    % Abstract superclass used by ns.EpochChannel and ns.TepochChannel
    % When working with a set of EpochChannel or TepochChannel objects one often wants to compute
    % derived measures (e.g. a spectrum from a signal), average in
    % different ways (across trials or channels, or subjects), or visualize
    % raw or computed data.
    %
    % This cache class prevents multiple round trips to the server to fetch
    % the data. Instead, the data are fetched once and stored internally in
    % a Matlab table (.T) that also tracks the various primary keys of each
    % row in the DJ table.
    %
    %  The instructions for a computation are provided as a struct
    %  that combines the function (fft, pspectrum,etc) and the (optional) 
    %  input arguments. Each of the fields of the struct specifies one
    %  operation; they are executed in order. 
    %
    % SEE ns.cachce/compute for examples
    %
    % EXAMPLES
    % Determine the power spectral density between 0.5 and 50 Hz of a set of epochs 
    %  e=ns.EpochChannel 
    %  plot(e)  ; % Plots epoch time courses averaged acros trials & channels.
    % Compute the power spectrum with pspectrum
    %  fun.pspectrum = struct('FrequencyLimits',[0.5 50])
    %  G = compute(e,fun)
    % Plot the results (per condition)   
    %  clf;for g= 1:height(G),plot(G.frequency{g},G.power{g});hold on;end;legend(G.condition)
    %
    %
    % BK - Dec 2025
    % MO - 2026
    properties (Constant)
        GROUPVARS = ["subject" "session_date" "starttime" "paradigm" "condition" "trial" "channel"];
    end
    properties (GetAccess =public,SetAccess = protected)
        T (:,:) table  = table;  % The Matlab table that stores the data
        qry (1,1) string =""     % The query that fetched the data
        independent (1,:) string = "time" % The name(s) of the independent variables
        dependent (1,:) string = "signal"  % The name(s) of the dependent variables
        time (1,:) double = [] % Computed by fill
        samplingRate (1,1) double =NaN % Computed by fill
    end

    methods
        function v =get.T(o)
            % Fill the cache and return as as table
            fill(o);
            if ismember("dependent", o.T.Properties.VariableNames)
                % Tepoch - rename
                o.dependent = o.T.dependent(1);
                o.independent = o.T.independent(1);
                o.T = renamevars(o.T,["signal" "x"], [o.dependent o.independent]);
                o.T = removevars(o.T,["dependent" "independent"]);
            end
            v = o.T;
        end
    end

    methods (Access = public)
        function plot(o,pv)
            % Plot y as a function of x for all rows in the table.
            % Set the 'average' input to select which aspects to average
            % over. By default, trials and channels are averaged.
            arguments
                o (1,1) {mustHaveRows}
                pv.delta (1,1) string = ""              % Show the difference using this named condition as the reference
                pv.channel (:,1) double = []            % Select a subset of channels
                pv.trial (:,1) double = []            % Select a subset of trials
                pv.average (1,:) string {mustBeMember(pv.average,["starttime" "condition" "trial" "channel" "subject" "session_date" ""])} = ["trial" "channel"]  % Average over these dimensions
                pv.tilesPerPage (1,1) double = 6        % Select how many tiles per page.
                pv.linkAxes (1,1) logical = false        % Force the same xy axes on all tiles in a figure
                pv.raster (1,:) string = ""            % Set to true to show trials as rasters (removes "trial" from pv.average)
                pv.newTileEach = ["paradigm" "subject" "session_date" "starttime"];  % Start a new tile when any of these parameters change.
                pv.figure = []  % Creates new figures if empty.
            end

            %% Fill the cache, then perform averaging per group
            fill(o);% Fill the cache if needed
            dimension = unique(o.T.dimension);
            
           % Raster plot cannot average over 
           pv.average = setdiff(pv.average,pv.raster,'stable');
           
            grouping = setdiff(ns.cache.GROUPVARS,[pv.average pv.raster],'stable');

            % Epochs always contain signal and time
            xName = o.independent;
            yName = o.dependent;
            G = compute(o,struct("msten",[]),x=xName,y=yName,average= pv.average,channel=pv.channel,trial=pv.trial);            
            x = G{1,xName};
            if xName =="time" && numel(x) ==3
                x = linspace(x(1),x(2),x(3))';
            end
            if pv.raster ~=""
                % Concatenate the trials into a raster matrix in G.
                rasterGrouping = setdiff(ns.cache.GROUPVARS,[pv.raster pv.average]);
                P = groupsummary(G, rasterGrouping, @(x) x(1,:), ["align" xName "paradigm"]);
                P = renamevars(P,["fun1_align" "fun1_"+xName "fun1_paradigm"],["align" xName "paradigm"]);
                G = groupsummary(G,rasterGrouping,@(x) ({cat(1,x)}),["mean" "ste" "n"]);
                G = renamevars(G,["fun1_mean" "fun1_ste" "fun1_n"],["mean" "ste" "n"]);
                G = innerjoin(G,P);
                pv.newTileEach = union(pv.newTileEach,"condition");
                G= sortrows(G,intersect([ "subject" "session_date" "starttime" "condition" "channel" "trial" "paradigm"],G.Properties.VariableNames,'stable'));
            end


            %% Figure
            tileCntr=0;
            nrTimeSeries = height(G);
            pv.newTileEach = intersect(pv.newTileEach,G.Properties.VariableNames);
            legStr =string([]);
            for i = 1:nrTimeSeries
                align =G.align(i);
                if i==1 || (~isempty(pv.newTileEach) && any(G{i,pv.newTileEach} ~= G{i-1,pv.newTileEach}))
                    % New  subject, session or experiment in a new tile
                    if i>1 &&  ~isempty(legStr)
                        % Add the legend string to the existing tile
                        legend(h,legStr);
                    end
                    if isempty(pv.figure)
                    if mod(tileCntr,pv.tilesPerPage)==0 
                        if i>1 && pv.linkAxes
                            linkaxes(gcf().Children().Children())
                        end
                        figure;
                    end
                    else
                        figure(pv.figure); % Add to existing
                    end
                    % Start a new tile with empty handles
                    nexttile;
                    tileCntr =tileCntr+1;
                    h = [];
                    legStr = string([]);
                    hold on
                end
                if pv.raster~=""
                    % Show each condition in a separate tile
                    nrTrials= size(G.mean{i},1);
                    imagesc(x,1:nrTrials, cell2mat(G.mean{i}')')
                    axis xy
                    n = mean(cell2mat(G.n{i}'),"all");
                    ylabel ("")
                    titlePV= setdiff(["paradigm" rasterGrouping],"",'stable');
                    ttlStr = strjoin(string(G{i,titlePV}),"/");
                else
                    m = G.mean{i,:};
                    ste = G.ste{i,:};
                    n = mean(G.n{i,:});                    
                    h = [h plot(x,m)];                %#ok<AGROW>
                    p = patch([x ;  flip(x)]',[m+ste ; flip(m-ste)]',h(end).Color,FaceAlpha= 0.5);
                    p.EdgeColor = h(end).Color;
                    plot(xlim,[0 0],'k');
                    ylabel 'EP (\muV)'
                    legStr = [legStr dimension + "=" + G.condition(i)]; %#ok<AGROW>
                    titlePV= setdiff(["paradigm" grouping],"condition",'stable');
                    ttlStr = strjoin(string(G{i,titlePV}),"/");
                end
                title (ttlStr + " (n=" + string(n) +")",'Interpreter','none');
                xlabel (sprintf('Time after %s.%s (s)',align.plugin,align.event));
                % If delta is not empty, add the difference wave.
                if pv.delta ~="" && G.condition(i) ~=pv.delta
                    matchVars = setdiff(grouping,"condition");
                    if isempty(matchVars)
                        matchG = G;
                    else
                        matchG = innerjoin(G,G(i,matchVars));
                    end
                    reference = find(matchG.condition ==pv.delta);
                    if ~isempty(reference)
                        y = m - matchG.mean{reference,:};
                        ste = ste +matchG.ste{reference,:};
                        h = [h plot(x,y)];                     %#ok<AGROW>
                        p = patch([x;flip(x)]',[y+ste;flip(y-ste)]',h(end).Color,FaceAlpha= 0.5);
                        p.EdgeColor = h(end).Color;
                        legStr = [legStr G.condition(i)+"-"+ pv.delta]; %#ok<AGROW>
                    end
                end
            end
            if ~isempty(legStr); legend(h,legStr); end % For the last tile
            if pv.linkAxes
                warning('off','MATLAB:linkaxes:RequireDataAxes')
                linkaxes(gcf().Children().Children())
                warning('on','MATLAB:linkaxes:RequireDataAxes')
            end
        end

        function [G,dv,idv] =  compute(o,fun,pv)
            % Compute derived measures from the EpochChannel table.
            % fun - struct with fields corresponding to one of the 
            % functions listed below. FFT and pspectrum take a struct
            % containing named options; all spectral functions use this
            % struct-based option format.
            %   .fft:       Uses `fft()` to compute amplitude, phase as a function
            %                   of frequency. Supports the named `n` option.
            %                   Sampling rate is determined automatically.
            %                EXAMPLE: fun.fft = struct('n',128)
            %   .pspectrum: Uses `pspectrum()` to compute power spectral
            %                   density as a function of frequency
            %                   Sampling rate is determined automatically.
            %               EXAMPLE : fun.pspectrum = struct('FrequencyLimits',[0 50])
            %   .pmtm:      Uses the multitaper pmtm function. Options are
            %               supplied as a struct with these fields:
            %               tapertype: 'slepian' (default) or 'sine'
            %               nw:        4 (default) or a positive scalar
            %               m:         7 (default), an integer scalar, or vector
            %               nfft:      [] (default) or an integer
            %               fs:        sample rate, a positive scalar
            %               f:         frequencies, a vector
            %               EXAMPLE:
            %               fun.pmtm = struct('nw',4,'nfft',256,'fs',250);
            %   .msten:     Mean,standard error, and N (used by the ns.EpochChannel/plot function)
            %               EXAMPLE fun.msten = {};
            %   .wavelet:   Wavelet spectrogram using the FWHM approach.
            %               Sampling rate is supplied automatically. The
            %               options are supplied as a struct. Its fields
            %               specify the number of frequencies ('nfrex'),
            %               frequency limits ('limits'), and the FWHM at
            %               the lowest and highest frequencies ('fwhm').
            %               EXAMPLE:
            %               fun.wavelet = struct('nfrex',40,'fwhm',[2 0.2], ...
            %                   'limits',[0.5 50]);
            %
            %   .snr:       Calculate SNRs as a ratio between the power
            %                  at a given frequency and the average power
            %                  at its neighboring (noise) frequencies
            %                  excluding immediate neighbors. 
            %               Options are supplied as a struct with fields
            %               signalHalfWidth and noiseHalfWidth, both in Hz.
            %
            %               for a given frequency f_i, the noise power
            %                  is the average power in the range of
            %                   f_i + [1,-1].*noiseHalfwidth
            %                  excluding
            %                   f_i + [1,-1].*signalHalfwidth
            %               signalHalfWdith must be smaller than
            %               noiseHalfWidth and the signalHalfWidth must be
            %               bigger than the frequency resolution of the
            %               spectral analysis.
            %
            %               Note that this compute fun needs others to
            %               work properly, for instance:
            %                   fun.pspectrum = {'FrequencyLimits',[0 50]};
            %                   fun.snr = struct('signalHalfWidth',1, ...
            %                       'noiseHalfWidth',2);
            %               This will first compute the power spectrum and
            %               then the snr for all the frequencies in that
            %               spectrum. With pmtm this can be finetuned to only
            %               compute the power at specific frequencies of
            %               interest:
            %                   fun.pmtm = {4,1:24,250}; % compute power at 1:24 Hz
            %                   fun.snr = struct('signalHalfWidth',2, ...
            %                       'noiseHalfWidth',4);
            %               Or you can compute the snr at all frequencies
            %               and then find the peaks in the snr that are
            %               close to a set of frequencies of interest:
            %                   fun.pspectrum   = {'FrequncyLimits',[0 50]}; 
            %                   fun.snr = struct('signalHalfWidth',1, ...
            %                       'noiseHalfWidth',2); % Determine SNR
            %                   fun.peak = struct('searchFrequencies', ...
            %                       [2 6 10 24], 'searchRangeHalfWidth',1);
            %                   snr within 1 Hz from 2,6, 10,24 Hz.
            %            
            %             
            %   .peak:      Finds peak locations and magnitudes around
            %                   specific frequencies within a search window.
            %               Options are supplied as a struct with fields
            %               searchFrequencies and searchRangeHalfWidth.
            %
            %
            % channel  - Select a subset of channels
            % trial    - Select a subset of trials
            % timeWindow - Select a time window
            % average  - Analysis is applied after averaging over these
            %               fields.
            %               Defaults to ["trial" "channel"], but can be
            %               any of ["subject" "session_date" "starttime" "condition" "trial" "channel"]
            %           Use "" to compute for each trial and channel.
            % OUTPUT
            % G  - A table with the results
            % dv -  The name of the dependent variable.
            % idv - the name of the independent variable.
            arguments
                o (1,1)
                fun  (1,1) struct
                pv.channel (:,1) double = []            % Select a subset of channels
                pv.trial (:,1) double = []            % Select a subset of trials
                pv.timeWindow (1,2) double = [-inf inf]  % Select a time window to operate on
                pv.average (1,:) string {mustBeMember(pv.average,["" "subject" "session_date" "starttime" "condition" "trial" "channel"])} = ["trial" "channel"]
                pv.x (1,1) string = o.independent  % Name of the independent variable
                pv.y (1,1) string = o.dependent    % Name of the dependent variable
            end
            fill(o);% Fill the cache
            idv = pv.x;
            dv = pv.y;

            %% Restrict the T by function input args and time window selection
            stay = true(height(o.T),1);
            if ~isempty(pv.channel)
                % Restrict channels for this operation
                stay = stay & ismember(o.T.channel,pv.channel);
            end
            if ~isempty(pv.trial)
                % Restrict trials for this operation
                stay = stay  & ismember(o.T.trial,pv.trial);
            end
            restrictedT = o.T(stay,:);
            if isempty(restrictedT)
                error('No data in this table');
            end
            if any(isfinite(pv.timeWindow))
                assert(ismember("time",o.T.Properties.VariableNames),"timeWindow restriction can only be used on a cache with a time column.")
                % Crop to the timeWindow for this operation.
                t = restrictedT{1,"time"}; % Time in secods (Taking first row as all should be the same)
                if numel(t)==3
                    t = linspace(t(1),t(2),t(3));
                end
                keep = do.ifwithin(t,pv.timeWindow/1000);
                restrictedT(:,dv)= rowfun(@(x) {x(keep,:)},restrictedT(:,dv) ,'ExtractCellContents',true);
                t= t(keep);
                assert(~isempty(t),'No time points left in the analysis window ([%f %f])',pv.timeWindow(1),pv.timeWindow(2));
                restrictedT.time = repmat([t(1) t(end) numel(t)],height(restrictedT),1);
            end

            %%  Average/group
            if ~(isscalar(pv.average) && pv.average=="")
                if ismember("condition",pv.average) && ~ismember("trial",pv.average)
                    pv.average = [pv.average "trial"];
                end
                grouping = setdiff(ns.cache.GROUPVARS,pv.average,'stable');
                
                [grp,G] = findgroups(restrictedT(:,grouping));
                % Average per group
                if isfield(fun,"msten")
                    % Special case; caller asks for the mean only (and ste
                    % and n)
                    M = splitapply(@ns.cache.do_msten,restrictedT.(dv),grp);                    
                else
                    % Average signal then that will be processed by the fun
                    % below.
                    if iscell(restrictedT{1,dv})
                        M = splitapply(@(x) {mean(cat(2,x{:}),2,"omitmissing")},restrictedT.(dv),grp);
                    else
                        M = splitapply(@(x) {mean(x,1,"omitmissing")'},restrictedT.(dv),grp);
                    end
                end

                % Combine with align/time/paradigm information. Note this
                % assumes these are constant across the group (picking
                % only the first here). fill() assures this is the case.
                P = groupsummary(restrictedT, grouping, @(x) x(1,:), ["align" idv]);                
                P = renamevars(P,["fun1_align" "fun1_"+idv ],["align" idv ]);

                % add trial counts
                nT = groupsummary(restrictedT, grouping, @(x) numel(unique(x)), "trial");
                nT = renamevars(nT,"fun1_trial", "nrtrials");

                % add channel counts
                nCh = groupsummary(restrictedT, grouping, @(x) numel(unique(x)), "channel");
                nCh = renamevars(nCh,"fun1_channel", "nrchannels");

                G = innerjoin(G,P); 
                G = innerjoin(G,nT); 
                G = innerjoin(G,nCh); 
                G = removevars(G, "GroupCount");
            else
                % No averaging. Just put the signal into M
                G = restrictedT;
                G.nrtrials = ones(height(G),1);
                G.nrchannels = ones(height(G),1);
                M = restrictedT.(dv);
            end
            nrGrps = height(M);

            %% Determine which function to compute
            % Map string to function handle and do error checking
            if isfield(fun,"msten")
                % it outputs multiple dv
                dv = ["mean", "ste", "n"];
                idv = "time";
                % Combine with G
                G = [G M];
            else
                % Compute one or more functions
                funs = fieldnames(fun);
                nrFuns = numel(funs);
                fprintf('Applying %d functions (%s) to %d elements\n',nrFuns,strjoin(funs,"/"),size(M,1))
                for f = 1:nrFuns
                    thisFun = funs{f};
                    if isstruct(fun.(thisFun))
                        thisOptions = namedargs2cell(fun.(thisFun));
                    else
                        thisOptions = {};
                    end
                    data = {M}; % input to the function called by splitapply
                    srate= o.samplingRate; % Local copy to avoid broadcasting o in the parfor
                    switch thisFun
                        case "fft"
                            funN = @(x) ns.cache.do_fft(x,srate,thisOptions{:});
                            dv = ["amplitude" "phase"];
                            idv = "frequency";
                        case "pspectrum"
                            funN = @(x) ns.cache.do_pspectrum(x,srate,thisOptions{:});
                            idv = "frequency";
                            dv = "power";
                        case "pmtm"
                            funN = @(x) ns.cache.do_pmtm(x,'fs',srate,thisOptions{:});
                            idv = "frequency";
                            dv = "power";
                        case "wavelet"                         
                            funN = @(x) ns.cache.do_wavelet(x,srate,thisOptions{:});
                            % do_wavelet() does not output time
                            idv = "frequency";
                            dv = "power";                        
                        %% Cases below take G (the result of previous computation) as their input
                        case "snr"
                            assert(nrFuns>1,"snr cannot run on its own; fun needs a spectral power estimate.");
                            funN = @(varargin) ns.cache.do_snr(varargin{:},thisOptions{:});
                            data{end+1} = G.(idv); %#ok<AGROW> % The IDV of the previous comp is passed to the function; this should be a set of frequencies.
                            idv = "frequency";
                            dv = "snr";                        
                        case 'peak'
                            assert(nrFuns>1,"peak cannot run on its own; fun needs a spectral power or snr estimate.");
                            funN = @(varargin) ns.cache.do_search_peaks(varargin{:}, thisOptions{:});
                            data{end+1} = G.(idv); %#ok<AGROW> % frequency
                            idv = "searchFrequency";
                            dv = ["peakFrequency","magnitude"];
                        otherwise
                            error('Unknown function %s', thisFun);
                    end

                    %% Apply the fun to the mean signal
                    pool = nsParPool;
                    if isempty(pool)
                        x = splitapply(funN,data{:},(1:nrGrps)');
                    else
                        xCell = cell(nrGrps,1);
                        progressQueue = parallel.pool.DataQueue;
                        progressCount = 0;
                        progressListener = afterEach(progressQueue,@updateProgress); %#ok<NASGU>
                        fprintf('Applying %s in parallel: 0/%d',thisFun,nrGrps);
                        parfor iGrp = 1:nrGrps                            
                            groupData = cellfun(@(d) d(iGrp,:),data,UniformOutput=false);
                            xCell{iGrp} = funN(groupData{:}); %#ok<PFBNS>
                            send(progressQueue,1);
                        end
                        x = vertcat(xCell{:});
                        fprintf('\n');
                    end
                    % Combine with G, if G and x have common variable
                    % names, x overwrites G
                    G = ns.cache.horzcat_results_(G, x);

                    % change M for the next computation
                    M = table2cell(G(:,dv));
                end
            end
            % Sort in consistent order - not matched to the tbl query
            G= sortrows(G,intersect(["subject" "session_date" "starttime" "paradigm"  "condition" "channel" "trial"],G.Properties.VariableNames,'stable'));

            function updateProgress(~)
                progressCount = progressCount + 1;
                progressStep = max(1,ceil(nrGrps/10));
                if mod(progressCount,progressStep) == 0 || progressCount == nrGrps
                    fprintf('\rApplying %s in parallel: %d/%d',thisFun,progressCount,nrGrps);
                end
            end
        end
    end



    methods (Static)
        % Compute functions that take a signal with some options and return
        % a table with one or more output columns. Note that each column
        % should contain a row vector of results.
        function v = do_fft(signal, fs, pv)
        % do_fft - Computes FFT amplitude and phase for each
        %               epoch. Only includes real frequencies.
        %
        % Outputs (table columns):
        %   amplitude: Amplitude of the FFT.
        %   phase: Phase of the FFT.
        %   frequency: Corresponding real frequencies.
            arguments
                signal cell
                fs (1,1) double {mustBeFinite,mustBePositive}
                pv.n double {mustBePositiveIntegerOrEmpty} = []
            end        

            signal =cat(2,signal{:}); % Concatenate epochs
            % Compute FFT for each slice along time dim 1
            if ~isempty(pv.n)
                fftResult = fft(signal,pv.n);
            else
                fftResult = fft(signal);
            end

            % Calculate real frequencies
            N = size(fftResult, 1);
            if mod(N, 2) == 0
                freq = (0:N/2) * fs / N;
                idx = 1:N/2+1;
            else
                freq = (0:(N-1)/2) * fs / N;
                idx = 1:(N+1)/2;
            end

            amplitude = 2*abs(fftResult(idx,:,:)/sqrt(N));
            phase = angle(fftResult(idx,:,:));
            % Return as table with results as row vectors
            v= table(amplitude',phase',freq,'VariableNames',{'amplitude','phase','frequency'});
        end
        function v = do_pspectrum(signal, fs, pv)
            % Compute power spectral density using MATLAB's pspectrum function.
            arguments
                signal cell
                fs (1,1) double {mustBeFinite,mustBePositive}
                pv.FrequencyLimits (1,2) double {mustBeFrequencyLimitsOrEmpty} = [0 fs/2]
                pv.FrequencyResolution double {mustBePositiveScalarOrEmpty} = []
                pv.Leakage double {mustBeUnitIntervalOrEmpty} = [0.5]
                pv.MinThreshold double {mustBeScalarOrEmpty} = []
                pv.Reassign logical {mustBeScalarOrEmpty} = []
                pv.TwoSided logical {mustBeScalarOrEmpty} = []
            end
          
            % Table with power and frequency
            signal =cat(2,signal{:}); % Concatenate epochs
            signal = signal - mean(signal,1,"omitmissing");
            options = namedargs2cell(pv);
            for iOption = numel(options)-1:-2:1
                if isempty(options{iOption+1})
                    options(iOption:iOption+1) = [];
                end
            end
            [power, freq] = pspectrum(signal, fs, 'power', options{:});
            v= table(power',freq','VariableNames',{'power','frequency'});
        end
        function v = do_pmtm(signal,pv)
            arguments
                signal cell
                pv.tapertype (1,1) string {mustBeMember(pv.tapertype,["slepian" "sine"])} = "slepian"
                pv.nw (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 4
                pv.m double {mustBeSineTaperOption} = 7
                pv.nfft double {mustBePositiveIntegerOrEmpty} = []
                pv.fs (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
                pv.f double {mustBeVectorOrEmpty} = []
            end
            % Multitaper power and frequency
            signal =cat(2,signal{:});
            signal(isinf(signal) | isnan(signal))=0;
            assert(isempty(pv.nfft) || isempty(pv.f), ...
                "Specify only one of nfft and f.");
            frequencyInput = pv.f;
            if isempty(frequencyInput)
                frequencyInput = pv.nfft;
            end
            if pv.tapertype == "sine"
                [power, freq] = pmtm(signal,pv.m,'Tapers','sine', ...
                    frequencyInput,pv.fs);
            else
                [power, freq] = pmtm(signal,pv.nw,frequencyInput,pv.fs);
            end
            if isrow(freq) % make sure 1st dim is always frequency
                freq = freq';
                power = power';
            end
            % Make table, force rows
            v = table(power',freq','VariableNames',{'power','frequency'});
        end
        function v = do_wavelet(signal,fs, pv)
            arguments
                signal cell
                fs (1,1) double
                pv.nfrex (1,1) double = 40
                pv.fwhm (1,2) double = [2 0.2]
                pv.limits (1,2) double =[0.5 50]
            end
            % Code adapted from Cohen M. X. (2019). A better way to
            % define and describe Morlet wavelets for time-frequency
            % analysis. NeuroImage, 199, 81-86.
            % https://doi.org/10.1016/j.neuroimage.2019.05.048
            signal =cat(2,signal{:}); % Concatenate epochs
            nrSamples= size(signal,1);
            time = (0:nrSamples-1)/fs;
            % time-frequency parameters
            freq  = linspace(pv.limits(1),pv.limits(2),pv.nfrex)';
            fwhm = linspace(pv.fwhm(1),pv.fwhm(2),pv.nfrex)'; % variable fwhm
            assert(all(fwhm.*freq>=1),"The FWHM is too small (should have more than one cycle per window)");

            % setup wavelet and convolution parameters
            wavet = (-5:1/fs:5)';
            halfw = floor(length(wavet)/2)+1;
            nConv = nrSamples + length(wavet) - 1;
            % initialize time-frequency matrix
            spectrogram = zeros(pv.nfrex,nrSamples);
            % spectrum of data - for convolution with wavelets
            dataX = fft(signal,nConv);
            % loop over frequencies
            for fi=1:length(freq)
                % create wavelet
                waveX = fft( exp(2*1i*pi*freq(fi)*wavet).*exp(-4*log(2)*wavet.^2/fwhm(fi).^2),nConv );
                waveX = waveX./max(waveX); % normalize
                % convolve
                as = ifft( waveX.*dataX );
                % trim to valid part
                spectrogram(fi,:) = as(halfw+(1:nrSamples))';
            end
            power = abs(spectrogram).^2;
            % Store power spectrogram and frequency
            v = table({power'},freq',time,'VariableNames',{'power','frequency','time'});
        end
        function v = do_msten(x)
            % Mean, standard error, and N
            X =cat(2,x{:});
            v = {mean(X,2,"omitmissing")', ...  % Mean
                (std(X,0,2,"omitmissing")./sqrt(sum(~isnan(X),2,"omitmissing")))',...  % Standard error
                sum(~isnan(X),2,"omitmissing")'};  % Non-Nan N
            % Make a table.
            v = cell2table(v,"VariableNames",{'mean','ste','n'});
        end
        function v = do_snr(signal, freqs, pv)
            arguments
                signal (:,:)
                freqs (:,1)
                pv.signalHalfWidth (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 1
                pv.noiseHalfWidth (1,1) double {mustBeFinite,mustBeReal,mustBePositive}  = 2 
            end
            signalHalfWidth = pv.signalHalfWidth;
            noiseHalfWidth = pv.noiseHalfWidth;

            if iscell(signal)
                signal = cat(2,signal{:});
            end
            if iscell(freqs)
                freqs = cat(2,freqs{:});
            end
            assert(size(signal,1) == numel(freqs), "Signal and frequencies are of different length.");
            df = uniquetol(diff(freqs),1e-6); % frequency step
            assert(isscalar(df), "Frequencies are not regularly sampled.");
            % create the kernel
            halfWidth = floor(noiseHalfWidth/df);
            % must be even
            if rem(halfWidth,2), halfWidth = halfWidth + 1; end
            halfSkipWidth = floor(signalHalfWidth/df);
            if rem(halfSkipWidth,2), halfSkipWidth = halfSkipWidth + 1; end
            assert(signalHalfWidth<noiseHalfWidth,"Signal half width must be smaller than the noise half width");
            assert(signalHalfWidth>df,"Signal half width must be larger than the frequency spacing");
            assert(halfSkipWidth<halfWidth,"Signal and noise half widths do not leave any noise frequencies");
            kernel = ones(halfWidth,1);
            kernel(1:halfSkipWidth) = 0;
            kernel = [flip(kernel); 0; kernel];

            isFrequency0 = freqs == 0; % 0 Hz is only the DC offset
            signal(isFrequency0,:) = NaN; % DC offset should not be included in noise estimation
            % make signal log scale
            % noise is computed as mean instead of geomean that is more
            % appropriate for amplitudes. Log scaled signal mean acts
            % similar to geometric mean
            signal = log10(signal); 
            % Todo: THIS USED TO HAVE nan AS THIRD INPUT. wHY?
            noise = do.ndconv(signal, kernel)/sum(kernel); % conv is sum, make it mean
            snr = 10.^(signal - noise); % in log scale division becomes subtraction

            v = table(snr', freqs', VariableNames={'snr', 'frequency'});

        end
        function v = do_search_peaks(signal, freqs, pv)
            arguments
                signal (:,:)
                freqs (:,:)
                pv.searchFrequencies (1,:) double {mustBeFinite,mustBeReal,mustBeNonnegative}
                pv.searchRangeHalfWidth (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
            end
            searchFrequencies = pv.searchFrequencies;
            searchRangeHalfWidth = pv.searchRangeHalfWidth;

            if iscell(signal)
                signal = cat(2,signal{:});
            end
            if iscell(freqs)
                freqs = cat(2,freqs{:});
            end

            nSearchFreq = numel(searchFrequencies);
            nChannels = size(signal,2);
            [searchFreq, peakFreq, peakAmp] = deal(zeros(nSearchFreq,nChannels));
            for ii = 1:nSearchFreq
                
                sFreq = searchFrequencies(ii);
                searchWindow = sFreq + [-1, 1] .* searchRangeHalfWidth;
                isFrequency = do.ifwithin(freqs, searchWindow);
                frequencyIndices = find(isFrequency);
                if isempty(frequencyIndices)
                    error('No frequencies found in search window [%g, %g] Hz around %g Hz.', ...
                        searchWindow(1), searchWindow(2), sFreq);
                end
                [peakAmp(ii,:), peakIndices] = maxk(signal(isFrequency,:),1,1);

                peakFreq(ii,:) = freqs(frequencyIndices(peakIndices));
                searchFreq(ii,:) = sFreq;

            end

            v = table({searchFreq}, {peakFreq}, {peakAmp}, VariableNames={'searchFrequency', 'peakFrequency', 'magnitude'});
        end
    end

    methods (Static, Access = private)

        function result = horzcat_results_(G, M)
            % 1. Find overlapping names
            overlap = intersect(G.Properties.VariableNames, M.Properties.VariableNames);
            % 2. Remove them from G
            G = removevars(G, overlap);
            % 3. Horizontally concatenate
            result = [G, M];
        end
    end


    methods (Access= protected)
        function fill(o)
            % Fetch the data if the underlying query has changed
            [src] = getCacheQuery(o);
            if canonicalize(string(src.sql)) ~=canonicalize(o.qry)
                % Safety check; time and align should match for all rows
                % in the table.
                preFetch = fetchtable(src,'time','align');
                epochTime = preFetch.time;
                sameStart = isscalar(uniquetol(epochTime(:,1),1e-10));
                sameStop = isscalar(uniquetol(epochTime(:,2),1e-10));
                sameSamples = isscalar(unique(epochTime(:,3)));
                assert(sameStart && sameStop && sameSamples, ...
                    'Rows of the EpochChannel table must have identical start time, stop time, and number of samples.');
                assert(isscalar(unique({preFetch.align.plugin})),'Rows of the EpochChannel should be aligned to the same plugin.');
                assert(isscalar(unique({preFetch.align.event})),'Rows of the EpochChannel should be aligned to the same event.');
                o.T =fetchtable(src,'*','ORDER BY channel');
                o.qry = src.sql;
                o.time = linspace(epochTime(1,1),epochTime(1,2),epochTime(1,3));
                o.samplingRate  = epochTime(1,3)./(epochTime(1,2)-epochTime(1,1));            
            end

            function s = canonicalize(s)
                % 1. Find all aliases defined in AS clauses
                aliasPattern = '\s+AS\s+`?([$\w]+)`';
                s =  regexprep(s,aliasPattern,"AS ALIAS"); % name of the alias does not matter
                s = regexprep(s, '\s+', ' '); % single whitespace
                s = strtrim(s);
            end
        end
    end

    methods (Abstract, Access = protected)
        [src] = getCacheQuery(o)
    end

end

function mustBeSineTaperOption(value)
if isscalar(value)
    validateattributes(value,{'numeric'},{'real','finite','integer','positive'});
else
    validateattributes(value,{'numeric'},{'vector','real','finite'});
end
end

function mustBeVectorOrEmpty(value)
if ~isempty(value)
    validateattributes(value,{'numeric'},{'vector','real','finite'});
end
end

function mustBePositiveIntegerOrEmpty(value)
if ~isempty(value)
    validateattributes(value,{'numeric'},{'scalar','real','finite','positive','integer'});
end
end

function mustBePositiveScalarOrEmpty(value)
if ~isempty(value)
    validateattributes(value,{'numeric'},{'scalar','real','finite','positive'});
end
end

function mustBeUnitIntervalOrEmpty(value)
if ~isempty(value)
    validateattributes(value,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
end
end

function mustBeScalarOrEmpty(value)
if ~isempty(value)
    validateattributes(value,{'numeric','logical'},{'scalar','real'});
end
end

function mustBeFrequencyLimitsOrEmpty(value)
if ~isempty(value)
    validateattributes(value,{'numeric'},{'vector','numel',2,'real','finite','nondecreasing'});
end
end
