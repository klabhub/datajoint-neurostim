classdef (Abstract) cache < handle
    % Abstract superclass used by ns.EpochChannel and ns.TepochChannel
    %
    % When working with a set of EpochChannel or TepochChannel objects one 
    % often wants to compute derived measures (e.g. a spectrum from a signal), 
    % average in different ways (across trials or channels, or subjects), 
    % or visualize raw or computed data.
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
    % SEE ns.cachce/compute for detailed instructions and examples
    %
    % EXAMPLE
    % Determine the power spectral density between 0.5 and 50 Hz of a set of epochs
    %  e=ns.EpochChannel
    %  plot(e)  ; % Plots epoch time courses averaged acros trials & channels.
    % Compute the power spectrum with pspectrum
    %  fun.pspectrum = struct('FrequencyLimits',[0.5 50])
    %  G = compute(e,fun)
    % Plot the results (per condition)
    %  clf;for g= 1:height(G)
    %       plot(G.frequency{g},G.power{g});
    %       hold on;
    %       end;
    % legend(G.condition)
    %    
    % BK - Dec 2025
    % MO - 2026
    properties (Constant)
        AVERAGEVARS = ["subject" "session_date" "starttime" "paradigm" "condition" "trial" "channel"];
    end
    properties (GetAccess =public,SetAccess = protected)
        T (:,:) table  = table;  % The Matlab table that caches the data
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
                % Tepoch with named dv and idv - rename for easy reference
                o.dependent = o.T.dependent(1);
                o.independent = o.T.independent(1);
                o.T = renamevars(o.T,["signal" "x"], [o.dependent o.independent]);
                o.T = removevars(o.T,["dependent" "independent"]);
            end
            if o.T.Properties.VariableTypes(o.T.Properties.VariableNames==o.dependent) =="double"
                o.T.(o.dependent) = num2cell(o.T.(o.dependent),2); 
            end
            if o.T.Properties.VariableTypes(o.T.Properties.VariableNames==o.independent) =="double"
                o.T.(o.independent) = num2cell(o.T.(o.independent),2);
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
                o (1,1) 
                pv.delta (1,1) string = ""              % Show the difference using this named condition as the reference
                pv.channel (:,1) double = []            % Select a subset of channels
                pv.trial (:,1) double = []            % Select a subset of trials
                pv.average (1,:) string {mustBeMemberOrEmpty(pv.average,["starttime" "condition" "trial" "channel" "subject" "session_date" ""])} = ["trial" "channel"]  % Average over these dimensions
                pv.robust (1,1) logical = false;       % Set to true to use median and iqr
                pv.outlier (1,1) double = inf           % Tthreshold factor for outlier removal
                pv.tilesPerPage (1,1) double = 6        % Select how many tiles per page.
                pv.linkAxes (1,1) logical = false        % Force the same xy axes on all tiles in a figure
                pv.raster (1,:) string = ""            % Set to true to show trials as rasters (removes "trial" from pv.average)
                pv.line (1,1) logical = false           % Show offset lines instead of raster
                pv.newTileEach = ["paradigm" "subject" "session_date" "starttime" ];  % Start a new tile when any of these parameters change.
                pv.figure = []  % Creates new figures if empty.
                pv.xlim  (1,:) double = []
                pv.clim (1,:) double  = []
            end

            %% Fill the cache, then perform averaging per group
            fill(o);% Fill the cache if needed
            dimension = unique(o.T.dimension);

            % Raster plot cannot average over the dimension that should be rastered
            pv.average = setdiff(pv.average,pv.raster,'stable');
            % Determine the grouping variables
            grouping = setdiff(ns.cache.AVERAGEVARS,[pv.average pv.raster],'stable');

            % Epochs always contain signal and time
            xName = o.independent;
            yName = o.dependent;
            if isempty(pv.average)
                fun = struct();
            else
                fun.average = struct(...
                    "average",pv.average, ...
                    "robust",pv.robust, ...
                    "outlier",pv.outlier);
            end
            G = compute(o,fun,x=xName,y=yName,average=string.empty, ...
                channel=pv.channel,trial=pv.trial);
            if isempty(pv.average)
                plotXName = xName;
                plotYName = yName;
            else
                plotXName = "average_" + xName;
                plotYName = "average_" + yName;
            end
            
 
            if pv.raster ~=""
                % Concatenate the trials into a raster matrix in G.
                rasterGrouping = setdiff(ns.cache.AVERAGEVARS,[pv.raster pv.average]);
                P = groupsummary(G, rasterGrouping, @(x) x(1), ["align" plotXName "paradigm"]);
                P = renamevars(P,["fun1_align" "fun1_"+plotXName ],["align" plotXName ]);
                if isempty(pv.average)
                    G = groupsummary(G,rasterGrouping,@(x) {cat(1,x{:})},yName);
                    G = renamevars(G,"fun1_" +yName,yName);
                    G.n =repmat({1},height(G),1);
                else
                    averageVariables = ["average_" + yName, ...
                        "average_" + yName + "_error", ...
                        "average_" + yName + "_n"];
                    G = groupsummary(G,rasterGrouping,@(x) {cat(1,x{:})},averageVariables);
                    G = renamevars(G,"fun1_" + averageVariables,[plotYName plotYName + "_error" "n"]);
                end
                G = innerjoin(G,P);
                %pv.newTileEach = union(pv.newTileEach,"condition");
                G= sortrows(G,intersect([ "subject" "session_date" "starttime" "condition" "channel" "trial" "paradigm"],G.Properties.VariableNames,'stable'));
            end


            %% Figure
            tileCntr=0;
            nrTimeSeries = height(G);
            pv.newTileEach = intersect(pv.newTileEach,G.Properties.VariableNames);
            legStr =string([]);
            for i = 1:nrTimeSeries
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
                        drawnow
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
                    x = G.(plotXName){1}';
                    if xName =="time" && numel(x) ==3
                        x = linspace(x(1),x(2),x(3))';
                    end                    
                    nrTrials= size(G.(plotYName){i},1);                    
                    if pv.line
                        y =G.(plotYName){i};
                        y = y./max(abs(y),[],"all") + repmat((1:size(y,1))',[1 size(y,2)]);
                        plot(x,y)
                    else
                        imagesc(x,1:nrTrials, G.(plotYName){i})
                        axis xy
                        ylabel (pv.raster)
                        if ~isempty(pv.clim)
                            clim(pv.clim);
                        end
                    end
                    n = mean(G.n{i},"all");
                    hold on 
                    plot([0 0],ylim,'k')
                    keep = intersect(setdiff(ns.cache.AVERAGEVARS,pv.newTileEach),G.Properties.VariableNames);
                    ttlStr = strjoin(string(G{i,keep}),"/");
                elseif isempty(pv.average)
                    x = G.(plotXName){1}';
                    if xName =="time" && numel(x) ==3
                        x = linspace(x(1),x(2),x(3))';
                    end
                    y = G{i,plotYName};
                    if iscell(y)
                        y = cat(1,y{:});
                    end                    
                   
                    try
                    plot(x,y);
                    catch
                    end
                    hold on;
                    xlabel (plotXName)
                    ylabel (plotYName)
                    titlePV= setdiff(["paradigm" grouping],"condition",'stable');
                    ttlStr = strjoin(string(G{i,titlePV}),"/");
                    if ismember("condition",G.Properties.VariableNames)
                        legStr = [legStr dimension + "=" + G.condition(i)]; %#ok<AGROW>
                    end
                    n =NaN;
                else 
                    % An average has been determined
                    x =G.(plotXName){i}';
                    if xName =="time" && numel(x) ==3
                        x = linspace(x(1),x(2),x(3))';
                    end
                    m = G.(plotYName){i}';
                    err = G.(plotYName + "_error"){i}';
                    n = mean(G.(plotYName + "_n"){i});
                    if all(isnan(m)) || all(isnan(err))
                        % No data for this group, skip
                        warning("No data for this group, skipping");
                        continue
                    end
                    h = [h plot(x,m)];                %#ok<AGROW>
                    out = isnan(m);
                    p = patch([x(~out) ;  flip(x(~out))]',[m(~out)+err(~out) ; flip(m(~out)-err(~out))]',h(end).Color);
                    p.EdgeColor = h(end).Color;                    
                    p.FaceAlpha = 0.5;
                    plot(xlim,[0 0],'k');
                    ylabel (o.dependent,'Interpreter','none');
                    if ismember("condition",G.Properties.VariableNames)
                        legStr = [legStr dimension + "=" + G.condition(i)]; %#ok<AGROW>
                    end
                    titlePV= setdiff(["paradigm" grouping],"condition",'stable');
                    ttlStr = strjoin(string(G{i,titlePV}),"/");
                end
                title (ttlStr + " (n=" + string(n) +")",'Interpreter','none');
                xlabel (o.independent,'Interpreter','none');
                if ~isempty(pv.xlim)
                    xlim(pv.xlim)
                end
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
                        y = m - matchG.(plotYName){reference}';
                        err = err + matchG.(plotYName + "_error"){reference}';
                        h = [h plot(x,y)];                     %#ok<AGROW>
                        p = patch([x;flip(x)]',[y+err;flip(y-err)]',h(end).Color,FaceAlpha= 0.5);
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

        function [G,D,uGroup] =  compute(o,fun,pv)
            % Compute derived measures from the (T)EpochChannel table.
            % fun - struct with fields corresponding to one of the
            % functions listed below. 
            % Sampling rate is determined automatically, from the database.
            %   .fft:       Uses `fft()` to compute amplitude, phase as a function
            %                   of frequency. Supports the named `n` option.            
            %                EXAMPLE: fun.fft = struct()
            %   .pspectrum: Uses `pspectrum()` to compute power spectral
            %                   density as a function of frequency            
            %               EXAMPLE : fun.pspectrum = struct('FrequencyLimits',[0 50])
            %   .pmtm:      Uses the multitaper pmtm function. Options are
            %               supplied as a struct with these fields:
            %               tapertype: 'slepian' (default) or 'sine'
            %               nw:        4 (default) or a positive scalar
            %               m:         7 (default), an integer scalar, or vector
            %               nfft:      [] (default) or an integer            
            %               f:         frequencies, a vector
            %               EXAMPLE:
            %               fun.pmtm = struct('nw',4,'nfft',256);
            %   .average:    Average the result of the preceding function,
            %               including an error estimate and N. For example:
            %                   fun.pspectrum = struct('FrequencyLimits',[0 40]);
            %                   fun.average = struct('average',"trial");
            %               computes one spectrum per trial and then averages
            %               the spectra across trials. Set robust=true to use
            %               the median and IQR instead of the mean and STE.
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
            %                   fun.pspectrum   = struct('FrequencyLimits',[0 50]);
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
            % You can extend this functionality by specifying your own function handle as a field of the fun struct. 
            % The function must accept the signal and sampling rate as its first two arguments, and return a table with row vectors of results. 
            % The function must also return a second output, which is the name of the column in the table that is the independent variable.
            % For instance
            %        fun.myfun = @(signal,srate) myfun(signal,srate,'start',1,'mode','bla');
            %        The function myfun must return a table with one row per group and one column per dependent variable.
            %        Note how the additional inputs to your function are specified in the function handle.
            %
            % channel  - Select a subset of channels. Defaults to []: all channels
            % trial    - Select a subset of trials. Defaults to []: all trials
            %
            % Channel or trial seelection can be donw with a vector of 
            % channel/ trial numbers , or a function that takes the cached data 
            % table (o.T) as input and returns a logical vector that
            % indicates per row whether it should be included or not. 
            % For instance: trial   = @(T) (T.trial<100), will include only
            % trial numbers up to 100 in the computation.
            %
            % timeWindow - Select a time window
            % average  - Analysis is applied after averaging over these
            %               fields.
            %               Defaults to ["trial" "channel"], but can be
            %                any of ["subject" "session_date" "starttime" "condition" "trial" "channel"]
            %               Set to string.empty to avoid averaging.
            % OUTPUT
            % G  - A table with the results
            % D -  A dictionary mapping independent variables to dependent
            % variables (both of which are columns in G)
            % uGroup = the list of unique group names (or "_" if no
            % averaging was performed).
            arguments
                o (1,1)
                fun  (1,1) struct
                pv.channel (:,1)  {mustBeA(pv.channel,["double" "function_handle"])} = [] % Select a subset of channels
                pv.trial (:,1)  {mustBeA(pv.trial,["double" "function_handle"])} = [] % Select a subset of trials
                pv.timeWindow (1,2) double = [-inf inf]  % Select a time window to operate on
                pv.average (1,:) string {mustBeMemberOrEmpty(pv.average,["subject" "session_date" "starttime" "paradigm" "condition" "trial" "channel"])} = ["trial" "channel"]
                pv.robust (1,1) logical = false   % Set to true to determine median as average                
                pv.outlier (1,1) double = inf     % Threshold to remove outliers before averaging 
                pv.x (1,1) string = o.independent  % Name of the independent variable
                pv.y (1,1) string = o.dependent    % Name of the dependent variable
            end
            fill(o);% Fill the cache
            idv = pv.x;
            dv = pv.y;  
            srate= o.samplingRate; % Local copy to avoid broadcasting o in the parfor

            %% Restrict the T by function input args and time window selection
            stay = true(height(o.T),1);
            if ~isempty(pv.channel)
                % Restrict channels for this operation
                if isa(pv.channel,'function_handle')
                    % Evaluate the function
                    stay = stay & pv.channel(o.T);
                else
                    stay = stay & ismember(o.T.channel,pv.channel);
                end
            end

            if ~isempty(pv.trial)                
                % Restrict trials for this operation
                if isa(pv.trial,'function_handle')
                    % Evaluate the function
                    stay = stay & pv.trial(o.T);
                else
                    stay = stay  & ismember(o.T.trial,pv.trial);
                end
            end
            restrictedT = o.T(stay,:);
            if isempty(restrictedT)
                error('No data in this table');
            end
            
            if any(isfinite(pv.timeWindow))
                t = restrictedT.time{1}; % Time in seconds (all should be the same)              
                if numel(t)==3
                    t = linspace(t(1),t(2),t(3));
                end
                assert(ismember("time",o.T.Properties.VariableNames),"timeWindow restriction can only be used on a cache with a time column.")
                % Crop to the timeWindow for this operation.
               
                keep = do.ifwithin(t,pv.timeWindow/1000);
                if restrictedT.Properties.VariableTypes(string(restrictedT.Properties.VariableNames)==dv)=="double"
                    restrictedT.(dv) = restrictedT.(dv)(:,keep);
                else 
                    % Cell
                    restrictedT.(dv) = cellfun(@(x) x(:,keep), restrictedT.(dv),UniformOutput=false);
                end
                t= t(keep);
                assert(~isempty(t),'No time points left in the analysis window ([%f %f])',pv.timeWindow(1),pv.timeWindow(2));
                restrictedT.time = repmat({[t(1) t(end) numel(t)]},height(restrictedT),1);
            end

            %RestrictedT is a table with each subject/session/experiment/trial/channel as a row

            %%  Average/group
            if isempty(pv.average) 
                 % No averaging. Just put the signal into M
                G = restrictedT;
                G.nrtrials = ones(height(G),1);
                G.nrchannels = ones(height(G),1);
                M = restrictedT.(dv);                
                uGroup = "_";
                G.name = repmat(uGroup,height(G),1);
            else
                [G,~,averageValues,~,~] = ns.cache.averageData(...
                    restrictedT, restrictedT.(dv), pv.average, srate, pv.robust, pv.outlier);
                grouping = setdiff(ns.cache.AVERAGEVARS,pv.average,'stable');
                keep = intersect([grouping "align" idv "nrtrials" "nrchannels"], ...
                    string(G.Properties.VariableNames), 'stable');
                G = G(:,keep);
                varies = varfun(@(x) numel(unique(x)) > 1, G(:,grouping),OutputFormat='uniform');
                varies(grouping=="trial") = false; % Tracked separately
                if any(varies)
                    uGroup  = string(G{:,grouping(varies)});
                    if size(uGroup,2)>1
                        uGroup = join(string(uGroup), "/", 2);
                    end
                else
                    uGroup = repmat("all",height(G),1); % Averaging reduced this to a single group.
                end
                % M already contains the averaged signal for each group.
                M = averageValues;
                G.name = uGroup;
            end
            nrGrps = height(M);

            %% Determine which function to compute
            % Map string to function handle and do error checking
            D = containers.Map;
            % Compute one or more functions
            funs = fieldnames(fun);
            nrFuns = numel(funs);
            fprintf('Applying %d functions (%s) to %d elements\n',nrFuns,strjoin(funs,"/"),size(M,1))
            pool = nsParPool;
            progressQueue = [];
            progressListener = [];
            if ~isempty(pool)
                progressQueue = parallel.pool.DataQueue;
                progressListener = afterEach(progressQueue,@(~) ns.cache.advanceParforProgress());
            end
                for f = 1:nrFuns
                    thisFun = string(funs{f});
                    if isstruct(fun.(thisFun))
                        thisOptions = namedargs2cell(fun.(thisFun));
                    elseif isempty(fun.(thisFun)) || isa(fun.(thisFun),'function_handle')
                        thisOptions = {};
                    else
                        error('Function options for %s must be a struct, empty, or a function handle', thisFun);
                    end
                    if thisFun == "average"
                        % Special case, needs access to G
                        if isKey(D,idv)
                            D = remove(D,idv); % This will be replaced by average_idv
                        end
                        [G,M,D,idv,dv] = ns.cache.averageResults(...
                            G,M,D,idv,dv,srate,fun.(thisFun));
                        nrGrps = height(G);
                        continue
                    end

                    if ismember(thisFun,["snr" "peak"])
                        assert(f > 1 && ~isempty(idv), ...
                            '%s requires a preceding spectral function.',thisFun);
                    end

                    data = {M}; % Inputs indexed by group row                  
                    switch thisFun
                        case "fft"
                            funN = @(signal,srate) ns.cache.do_fft(signal,srate,thisOptions{:});                           
                        case "pspectrum"
                            funN = @(signal,srate) ns.cache.do_pspectrum(signal,srate,thisOptions{:});                            
                        case "pmtm"
                            funN = @(signal,srate) ns.cache.do_pmtm(signal,'fs',srate,thisOptions{:});                            
                        case "wavelet"
                            funN = @(signal,srate) ns.cache.do_wavelet(signal,srate,thisOptions{:});                                                        
                            %% Cases below take G (the result of previous computation) as their input
                        case "snr"
                            assert(nrFuns>1,"snr cannot run on its own; fun needs a spectral power estimate.");
                            funN = @(signal,freqs,srate) ns.cache.do_snr(signal,freqs,srate,thisOptions{:});
                            data{end+1} = G.(idv); %#ok<AGROW> % The IDV of the previous comp is passed to the function; this should be a set of frequencies.                           
                        case 'peak'
                            assert(nrFuns>1,"peak cannot run on its own; fun needs a spectral power or snr estimate.");
                            funN = @(signal,freqs,srate) ns.cache.do_search_peaks(signal,freqs,srate, thisOptions{:});
                            data{end+1} = G.(idv); %#ok<AGROW> % frequency                        
                        otherwise
                            if isa(fun.(thisFun),'function_handle')
                               % Check that the function takes two inputs (data and sampling rate) 
                               % returns two outputs (the table with results and the name of the idv column)
                               % The nargout==-1 is there to handle an
                               % anonymous function that uses deal to
                               % generate two outputs.
                               funN = fun.(thisFun);
                               assert(nargin(funN)==2,"The compute function (%s) must take two inputs (data and sampling rate)",thisFun);      
                               assert(nargout(funN)==2 || nargout(funN)==-1,"The compute function (%s) must return two outputs (the table with results and the name of the idv column)",thisFun);
                            else
                                error('Unknown function %s', thisFun);
                            end
                    end
                  
                    %% Apply the fun to the mean signal
                    % Temp cell to store results
                    xCell = cell(nrGrps,1);
                    idv  = repmat("",1, nrGrps);
                    if isempty(pool)
                        for iGrp = 1:nrGrps
                            groupData = cellfun(@(d) selectDataRow(d,iGrp),data,UniformOutput=false);
                            [xCell{iGrp},idv(iGrp)] = funN(groupData{:},srate);
                        end
                    else
                        ns.cache.resetParforProgress(thisFun,nrGrps);
                        parfor iGrp = 1:nrGrps
                            groupData = cellfun(@(d) selectDataRow(d,iGrp),data,UniformOutput=false);
                            [xCell{iGrp},idv(iGrp)] = funN(groupData{:},srate); %#ok<PFBNS>
                            send(progressQueue,1);
                        end
                        ns.cache.finishParforProgress();
                    end
                    R = vertcat(xCell{:}); % Results table
                    idv = unique(idv);
                    assert(isscalar(idv), ...
                        'Function %s returned inconsistent independent-variable columns.',thisFun);
                    R = renamevars(R,R.Properties.VariableNames,thisFun + "_" + R.Properties.VariableNames );
                    idv = thisFun + "_" + idv ;
                    dv = setdiff(R.Properties.VariableNames,idv); % EVerything but the idv

                    %% Check the format of the output 
                    % The dv is usuallly a row vector (or a scalar), but it has to be placed inside a cell to allow for some functions 
                    % that return a matrix. The only thing that is not allowed is a column vector.
                    % This is mainly here to avoid errors in the subsequent plotting function, which expects a row vector or a matrix.
                    for col = [dv idv]
                        assert(iscell(R{:,col}), 'The variable (%s) must be a cell array for each group', col);
                        assert(all(cellfun(@(v) isrow(v) || isscalar(v) || (size(v,1)>1 && size(v,2)>1), R{:,col})), 'The variable (%s) must be numerical values inside a cell) for each group', col);
                    end
                    axisSizes = cellfun(@numel,R{1,idv(1)});
                    for col = dv
                        valueSize = size(R{1,col}{1});
                        if numel(axisSizes)==1
                            matchesAxes = ismember(axisSizes,valueSize);
                        else
                            matchesAxes = isequal(sort(axisSizes),sort(valueSize));
                        end
                        assert(matchesAxes,"IDV and DV must match in size")
                    end
                    
                    % Pass the dependent variables to the next computation
                    M = table2cell(R(:,dv));        
                    
                    % Combine the results with G and store idv->dv mapping                    
                    D(idv)  = dv;                     
                    G = [G R]; %#ok<AGROW>                    
                end
                if ~isempty(pool)
                    delete(progressListener);
                end

            % Sort in consistent order - not matched to the tbl query
            G= sortrows(G,intersect(["subject" "session_date" "starttime" "paradigm"  "condition" "channel" "trial"],G.Properties.VariableNames,'stable'));
          
        end
    end



    methods (Static)
        % Compute functions that take a signal with some options and return
        % a table with one or more output columns. Note that each column
        % should contain a row vector of results.
        function [G,M,D,idv,dv] = averageResults(G,M,D,idv,dv,srate,options)
            % Average the output of a preceding compute function over cache dimensions.
            assert(isstruct(options) && isfield(options,"average"), ...
                'fun.average requires an ''average'' dimension, for example struct(''average'',"trial").');
            idv = string(idv);
            dv = string(dv);
            average = string(options.average);
            idvValues = G.(idv);
            [G,~,averageValues,errorValues,nValues,groupInfo] = ns.cache.averageData(...
                G,M,average,srate,getOption(options,'robust',false),getOption(options,'outlier',inf));
            ns.cache.validateIndependentVariable(idvValues,groupInfo.groupNumber,idv);
            idvValues = idvValues(groupInfo.first);
            grouping = groupInfo.grouping;
            varies = varfun(@(x) numel(unique(x)) > 1, G(:,grouping), ...
                OutputFormat='uniform');
            varies(grouping=="trial") = false;
            if any(varies)
                uGroup = string(G{:,grouping(varies)});
                if size(uGroup,2)>1
                    uGroup = join(string(uGroup), "/", 2);
                end
            else
                uGroup = repmat("all",height(G),1);
            end
            G.name = uGroup;
            oldVariables = intersect([idv dv], string(G.Properties.VariableNames), 'stable');
            G = removevars(G, oldVariables);
            idvOut = "average_" + idv;
            averageOut = "average_" + dv;
            errorOut = averageOut + "_error";
            nOut = averageOut + "_n";

            R = table(idvValues);
            R.Properties.VariableNames = cellstr(idvOut);
            for j = 1:numel(dv)
                R.(averageOut(j)) = averageValues(:,j);
                R.(errorOut(j)) = errorValues(:,j);
                R.(nOut(j)) = nValues(:,j);
            end

            G = [G R];
            M = averageValues;
            idv = idvOut;
            dv = averageOut;
            D(idv) = [dv errorOut nOut];
        end

        function [G,M,averageValues,errorValues,nValues,groupInfo] = averageData(G,M,average,srate,robust,outlier)
            % Group rows and average one or more cell/numeric data columns.
            average = string(average);
            mustBeMemberOrEmpty(average, ns.cache.AVERAGEVARS);
            variables = string(G.Properties.VariableNames);
            assert(all(ismember(average, variables)), ...
                'Cannot average over a dimension (%s) that is not present in the result table.',average);

            available = intersect(ns.cache.AVERAGEVARS, variables, 'stable');
            grouping = setdiff(available, average, 'stable');
            [groupNumber, ~] = findgroups(G(:,grouping));
            first = splitapply(@(x) x(1), (1:height(G))', groupNumber);
            nrGroups = numel(first);
            groupInfo = struct('first',first,'groupNumber',groupNumber,'grouping',grouping);

            nrTrials = [];
            if ismember("trial", average) && ismember("trial", variables)
                nrTrials = splitapply(@(x) numel(unique(x)), G.trial, groupNumber);
            end
            nrChannels = [];
            if ismember("channel", average) && ismember("channel", variables)
                nrChannels = splitapply(@(x) numel(unique(x)), G.channel, groupNumber);
            end

            G = G(first,:);
            remove = intersect(average, variables, 'stable');
            G = removevars(G, remove);
            if ~isempty(nrTrials), G.nrtrials = nrTrials; end
            if ~isempty(nrChannels), G.nrchannels = nrChannels; end

            if ~iscell(M)
                M = num2cell(M, 2);
            end
            nrColumns = size(M,2);
            averageValues = cell(nrGroups,nrColumns);
            errorValues = cell(nrGroups,nrColumns);
            nValues = cell(nrGroups,nrColumns);
            groupRows = accumarray(groupNumber,(1:numel(groupNumber))',[],@(x){x});
            for j = 1:nrColumns
                column = M(:,j);
                canUseFastMean = ~robust && isinf(outlier) && ...
                    all(cellfun(@(x) isnumeric(x) && isrow(x),column));
                if canUseFastMean
                    sizes = cellfun(@numel,column);
                    canUseFastMean = all(sizes == sizes(1));
                end
                if canUseFastMean
                    values = vertcat(column{:});
                    averageValues(:,j) = splitapply(@(x) {mean(x,1,'omitmissing')},values,groupNumber);
                    errorValues(:,j) = splitapply(@(x) {std(x,0,1,'omitmissing') ./ ...
                        sqrt(sum(~isnan(x),1,'omitmissing'))},values,groupNumber);
                    nValues(:,j) = splitapply(@(x) {sum(~isnan(x),1,'omitmissing')},values,groupNumber);
                else
                    for i = 1:nrGroups
                        result = ns.cache.do_average(column(groupRows{i}),srate, ...
                            robust=robust,outlier=outlier);
                        averageValues{i,j} = result.average{1};
                        errorValues{i,j} = result.error{1};
                        nValues{i,j} = result.n{1};
                    end
                end
            end
            M = averageValues;
        end

        function validateIndependentVariable(values,groupNumber,name)
            % Verify that an IDV is constant within every averaging group.
            if ~iscell(values)
                values = num2cell(values,2);
            end
            for i = 1:max(groupNumber)
                rows = find(groupNumber == i);
                reference = values{rows(1)};
                for j = 2:numel(rows)
                    candidate = values{rows(j)};
                    sameSize = isequal(size(reference),size(candidate));
                    if isnumeric(reference) && isnumeric(candidate) && sameSize
                        scale = max(1,max(abs(reference),[],'all'));
                        same = all(abs(reference-candidate) <= 1e-10*scale,'all');
                    else
                        same = isequaln(reference,candidate);
                    end
                    if ~same
                        error('ns:cache:IndependentVariableMismatch', ...
                            'Independent variable %s differs within averaging group %d.',name,i);
                    end
                end
            end
        end

        function [v,idv] = do_fft(signal, fs, pv)
            % do_fft - Computes FFT amplitude and phase for each
            %               epoch. Only includes real frequencies.
            %
            % Outputs (table columns):
            %   amplitude: Amplitude of the FFT.
            %   phase: Phase of the FFT.
            %   frequency: Corresponding real frequencies.
            arguments
                signal (:,1) {mustBeNumeric}
                fs (1,1) double {mustBeFinite,mustBePositive}
                pv.n double {mustBePositiveIntegerOrEmpty} = []
            end

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
            v= table({amplitude'},{phase'},{freq}, ...
                'VariableNames',{'amplitude','phase','frequency'});
            idv = "frequency";
        end
        function [v,idv] = do_pspectrum(signal, fs, pv)
            % Compute power spectral density using MATLAB's pspectrum function.
            arguments
                signal (:,1) {mustBeNumeric}
                fs (1,1) double {mustBeFinite,mustBePositive}
                pv.FrequencyLimits (1,2) double {mustBeFrequencyLimitsOrEmpty} = [0 fs/2]
                pv.FrequencyResolution double {mustBePositiveScalarOrEmpty} = []
                pv.Leakage double {mustBeUnitIntervalOrEmpty} = [0.5]
                pv.MinThreshold double {mustBeScalarOrEmpty} = []
                pv.Reassign logical {mustBeScalarOrEmpty} = []
                pv.TwoSided logical {mustBeScalarOrEmpty} = []
            end

            % Table with power and frequency
            signal = signal - mean(signal,1,"omitmissing");
            options = namedargs2cell(pv);
            for iOption = numel(options)-1:-2:1
                if isempty(options{iOption+1})
                    options(iOption:iOption+1) = [];
                end
            end
            out = isnan(signal);
            if any(out)
                fprintf(2,"Setting %.1f%% of the signal that to zero (removing NaN)\n",100*mean(out));
                signal(out)=0;
            end
            [power, freq] = pspectrum(signal, fs, 'power', options{:});
            v= table({power'},{freq'},'VariableNames',{'power','frequency'});
            idv = "frequency";
        end
        function [v,idv] = do_pmtm(signal,pv)
            arguments
                signal (:,1) {mustBeNumeric}
                pv.tapertype (1,1) string {mustBeMember(pv.tapertype,["slepian" "sine"])} = "slepian"
                pv.nw (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 4
                pv.m double {mustBeSineTaperOption} = 7
                pv.nfft double {mustBePositiveIntegerOrEmpty} = []
                pv.fs (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
                pv.f double {mustBeVectorOrEmpty} = []
            end
            % Multitaper power and frequency
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
            v = table({power'},{freq'},'VariableNames',{'power','frequency'});
            idv = "frequency";
        end
        function [v,idv] = do_wavelet(signal,fs, pv)
            arguments
                signal (:,1) {mustBeNumeric}
                fs (1,1) double
                pv.nfrex (1,1) double = 40
                pv.fwhm (1,2) double = [2 0.2]
                pv.limits (1,2) double =[0.5 50]
            end
            % Code adapted from Cohen M. X. (2019). A better way to
            % define and describe Morlet wavelets for time-frequency
            % analysis. NeuroImage, 199, 81-86.
            % https://doi.org/10.1016/j.neuroimage.2019.05.048
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
            v = table({power'},{freq',time},'VariableNames',{'power','xt'});
            idv = "xt";
        end
        function [v,idv] = do_average(signal,fs,pv)
            arguments
                signal
                fs (1,1) double %#ok<INUSA>                                
                pv.robust (1,1) logical = false                
                pv.outlier (1,1) double = inf
            end
            % Mean, standard error, and N
            if iscell(signal)
                signal =cat(1,signal{:});
            end
            if isfinite(pv.outlier)
                signal = rmoutliers(signal,"median","ThresholdFactor",pv.outlier);
            end
            if pv.robust
                av = median(signal,1,"omitmissing");
                err = arrayfun(@(j) iqr(signal(~isnan(signal(:,j)), j))/sqrt(sum(~isnan(signal(:,j)))),1:size(signal,2));                
            else
                av = mean(signal,1,"omitmissing");
                err= std(signal,0,1,"omitmissing")./sqrt(sum(~isnan(signal),1,"omitmissing"));
            end
            n = sum(~isnan(signal),1,"omitmissing");  % Non-Nan N
            % Make a table.
            v = table({av},{err},{n}, {1:size(av,2)}, ...
                'VariableNames',{'average','error','n','time'});
            idv = "time";  % time is not in v, but supplemented in the compute() code.
        end


        function [v,idv] = do_snr(signal, freqs, srate,pv)
            arguments
                signal (:,1)
                freqs (:,1)
                srate (1,1) double  %#ok<INUSA> %Not used but needed to match with other computes
                pv.signalHalfWidth (1,1) double {mustBeFinite,mustBeReal,mustBePositive} = 1
                pv.noiseHalfWidth (1,1) double {mustBeFinite,mustBeReal,mustBePositive}  = 2
            end
            signalHalfWidth = pv.signalHalfWidth;
            freqs = freqs(:);
            noiseHalfWidth = pv.noiseHalfWidth;

            assert(size(signal,1) == numel(freqs), "Signal and frequencies are of different length.");
            df = uniquetol(diff(freqs),1e-6); % frequency step
            assert(isscalar(df), "Frequencies are not regularly sampled.");
            % create the kernel to compute the noise
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
            noise = do.ndconv(signal, kernel,FillValue=NaN)/sum(kernel); % conv is sum, make it mean
            snr = 10.^(signal - noise); % in log scale division becomes subtraction

            v = table({snr'}, {freqs'}, VariableNames={'snr', 'frequency'});
            idv = 'frequency';
        end
        function [v,idv] = do_search_peaks(signal, freqs,srate, pv)
            arguments
                signal (:,1)
                freqs (:,1)
                srate (1,1) double  %#ok<INUSA> %Not used but needed to match with other computes           
                pv.searchFrequencies (1,:) double {mustBeFinite,mustBeReal,mustBeNonnegative}
                pv.searchRangeHalfWidth (1,1) double {mustBeFinite,mustBeReal,mustBePositive}
            end            
            freqs = freqs(:);            
            nSearchFreq = numel(pv.searchFrequencies);            
            [searchFreq, peakFreq, peakAmp] = deal(zeros(1,nSearchFreq));
            for ii = 1:nSearchFreq
                sFreq = pv.searchFrequencies(ii);
                searchWindow = sFreq + [-1, 1] .* pv.searchRangeHalfWidth;
                isFrequency = do.ifwithin(freqs, searchWindow);
                frequencyIndices = find(isFrequency);
                if isempty(frequencyIndices)
                    error('No frequencies found in search window [%g, %g] Hz around %g Hz.', ...
                        searchWindow(1), searchWindow(2), sFreq);
                end
                [peakAmp(ii), peakIndices] = maxk(signal(isFrequency,:),1,1);
                peakFreq(ii) = freqs(frequencyIndices(peakIndices));
                searchFreq(ii) = sFreq;
            end

            v = table({searchFreq}, {peakFreq}, {peakAmp}, ...
                VariableNames={'searchFrequency', 'frequency', 'magnitude'});
            idv = 'searchFrequency';
        end

        function resetParforProgress(label,total)
            state = struct(...
                'label', string(label), ...
                'total', total, ...
                'count', 0, ...
                'digits', strlength(string(total)), ...
                'messageLength', 0);
            setappdata(0,'nsCacheParforProgress',state);
            ns.cache.printParforProgress();
        end

        function advanceParforProgress()
            if ~isappdata(0,'nsCacheParforProgress')
                return
            end

            state = getappdata(0,'nsCacheParforProgress');
            state.count = state.count + 1;
            setappdata(0,'nsCacheParforProgress',state);
            ns.cache.printParforProgress();
        end

        function finishParforProgress()
            if ~isappdata(0,'nsCacheParforProgress')
                return
            end

            state = getappdata(0,'nsCacheParforProgress');
            state.count = state.total;
            setappdata(0,'nsCacheParforProgress',state);
            ns.cache.printParforProgress();
            fprintf('\n')
            rmappdata(0,'nsCacheParforProgress');
        end

        function printParforProgress()
            if ~isappdata(0,'nsCacheParforProgress')
                return
            end

            state = getappdata(0,'nsCacheParforProgress');
            backspace = repmat(sprintf('\b'),1,state.messageLength);
            message = sprintf('Applying %s in parallel: %*d/%d', ...
                state.label,state.digits,state.count,state.total);
            fprintf('%s%s',backspace,message);
            state.messageLength = strlength(message);
            setappdata(0,'nsCacheParforProgress',state);
            drawnow limitrate
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
                sameStart = isscalar(uniquetol(epochTime(:,1),0.1));
                sameStop = isscalar(uniquetol(epochTime(:,2),0.1));
                sameSamples = isscalar(unique(epochTime(:,3)));
                assert(sameStart && sameStop && sameSamples, ...
                    'Rows of the EpochChannel table must have identical start time, stop time, and number of samples.');
                assert(isscalar(unique({preFetch.align.plugin})),'Rows of the EpochChannel should be aligned to the same plugin.');
                assert(isscalar(unique({preFetch.align.event})),'Rows of the EpochChannel should be aligned to the same event.');
                o.T =fetchtable(src,'*','ORDER BY channel');              
                o.qry = src.sql;
                if ismember("signal",o.T.Properties.VariableNames) && iscell(o.T.signal) && all(cellfun(@iscolumn,o.T.signal))
                    o.T.signal = cellfun(@(x) (x'),o.T.signal,'UniformOutput',false);
                end
                o.time = linspace(epochTime(1,1),epochTime(1,2),epochTime(1,3));
                % Epoch endpoints are seconds; N samples span N-1 intervals.
                o.samplingRate = (epochTime(1,3)-1)/(epochTime(1,2)-epochTime(1,1));               
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

function value = getOption(options, name, default)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = default;
end
end

function d = selectDataRow(d,iGrp)
if iscell(d)
    d = d{iGrp};
    while iscell(d) && isscalar(d), d = d{1}; end
    if iscell(d), d = cat(2,d{:}); end
    d = d(:);
else
    d = d(iGrp,:);
    if isvector(d), d = d(:); end
end
end

function mustBeMemberOrEmpty(value, validValues)
if ~isempty(value)
    mustBeMember(value, validValues);
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
