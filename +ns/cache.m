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
                o (1,1) {mustHaveRows}
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
                pv.newTileEach = ["paradigm" "subject" "session_date" "starttime"];  % Start a new tile when any of these parameters change.
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
            G = compute(o,struct("averageforplot",''),outlier=pv.outlier, robust=pv.robust,x=xName,y=yName,average= pv.average,channel=pv.channel,trial=pv.trial);               
            
 
            if pv.raster ~=""
                % Concatenate the trials into a raster matrix in G.
                rasterGrouping = setdiff(ns.cache.AVERAGEVARS,[pv.raster pv.average]);
                P = groupsummary(G, rasterGrouping, @(x) x(1), ["align" xName "paradigm"]);
                P = renamevars(P,["fun1_align" "fun1_"+xName ],["align" xName ]);        
                if isempty(pv.average)
                    G = groupsummary(G,rasterGrouping,@(x) {cat(1,x{:})},"signal");
                    G = renamevars(G,"fun1_signal","signal");
                    G.n =repmat({1},height(G),1);
                else
                    G = groupsummary(G,rasterGrouping,@(x) {cat(1,x{:})},["average" "error" "n"]);
                    G = renamevars(G,["fun1_average" "fun1_error" "fun1_n"],["average" "error" "n"]);
                end
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
                    x = G.(xName){1}';
                    if xName =="time" && numel(x) ==3
                        x = linspace(x(1),x(2),x(3))';
                    end
                    if isempty(pv.average)
                        src = "signal";
                    else
                        src = "average";
                    end
                    nrTrials= size(G.(src){i},1);                    
                    if pv.line
                        y =G.(src){i};
                        y = y./max(abs(y),[],"all") + repmat((1:size(y,1))',[1 size(y,2)]);
                        plot(x,y)
                    else
                        imagesc(x,1:nrTrials, G.(src){i})
                        axis xy
                        ylabel (pv.raster)
                        if ~isempty(pv.clim)
                            clim(pv.clim);
                        end
                    end
                    n = mean(G.n{i},"all");
                    hold on 
                    plot([0 0],ylim,'k')
                    titlePV= setdiff(["paradigm" rasterGrouping],"",'stable');
                    ttlStr = strjoin(string(G{i,titlePV}),"/");
                elseif isempty(pv.average)
                    x = G.(xName){1}';
                    if xName =="time" && numel(x) ==3
                        x = linspace(x(1),x(2),x(3))';
                    end
                    y = G{i,yName};
                    if iscell(y)
                        y = cat(1,y{:});
                    end                    
                   
                    try
                    plot(x,y,'x');
                    catch
                    end
                    hold on;
                    xlabel (xName)
                    ylabel (yName)
                    ttlStr = "";
                    n =1;
                else 
                    % An average has been determined
                    x =G.(xName){i}';
                    if xName =="time" && numel(x) ==3
                        x = linspace(x(1),x(2),x(3))';
                    end
                    m = G.average{i}';
                    err = G.error{i}';
                    n = mean(G.n{i});
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
                        y = m - matchG.signal{reference,:};
                        err = err +matchG.error{reference,:};
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
            %   .average:     Average, an error estimate, and N (used by the ns.EpochChannel/plot function)
            %               EXAMPLE fun.average= {}; 
            %                           Returns mean and ste  
            %                      fun.average = {'robust',true}           
            %                           returns the median and iqr
            %               EXAMPLE fun.robust_average = {};
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
            % channel  - Select a subset of channels
            % trial    - Select a subset of trials
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
                pv.channel (:,1) double = []            % Select a subset of channels
                pv.trial (:,1) double = []            % Select a subset of trials
                pv.timeWindow (1,2) double = [-inf inf]  % Select a time window to operate on
                pv.average (1,:) string {mustBeMemberOrEmpty(pv.average,["subject" "session_date" "starttime" "condition" "trial" "channel"])} = ["trial" "channel"]
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
                G.group = repmat(uGroup,height(G),1);
            else                
                %if ismember("condition",pv.average) && ~ismember("trial",pv.average)
                %    pv.average = [pv.average "trial"];
                %end
                grouping = setdiff(ns.cache.AVERAGEVARS,pv.average,'stable');

                [grp,G] = findgroups(restrictedT(:,grouping));
                varies = varfun(@(x) numel(unique(x)) > 1, G(:,grouping),OutputFormat='uniform');     
                varies(grouping=="trial") = false; % Tracked separately
                if any(varies)
                    uGroup  = G{:,grouping(varies)};
                    if size(uGroup,2)>1
                        uGroup = join(string(uGroup), "/", 2);
                    end
                else
                    uGroup = repmat("all",height(G),1); % Averaging reduced this to a single group.
                end
                % Average signal per group
                if isfield(fun,"averageforplot")
                    % Special case; averaging for the plot function
                    % (includes a variance estimate)
                   M = splitapply(@(x) ns.cache.do_average(x,srate,robust=pv.robust,outlier=pv.outlier),restrictedT.(dv),grp);
                else
                    % Average for the fun computed below.
                    if pv.robust
                        av = @median;
                    else
                        av = @mean;
                    end
                    if isfinite(pv.outlier)
                        % Average after outlier removal
                        M = splitapply(@(x) {av(rmoutliers(cat(1,x{:}),"median","ThresholdFactor",pv.outlier),1,"omitmissing")}, ...
                            restrictedT.(dv),grp);
                    else
                        % Average
                        M = splitapply(@(x) {av(cat(1,x{:}),1,"omitmissing")}, ...
                            restrictedT.(dv),grp);
                    end
                end                                
                
                % Combine with align/time (or other idv) information. Note this
                % assumes these are constant per group (picking
                % only the first here). fill() assures this is the case.
                P = groupsummary(restrictedT, grouping, @(x) (x(1,:)), ["align" idv]);
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
                G.group = uGroup;              
            end
            nrGrps = height(M);

            %% Determine which function to compute
            % Map string to function handle and do error checking
            D = containers.Map;
           if isfield(fun,"averageforplot")
                % The average has already been determined above; just combine with G
                if iscell(M)
                    % No average
                    G.signal = M;
                else
                    % Averaged
                    G = [G M(:,setdiff(M.Properties.VariableNames,"time"))];
                end
                D('time') = {'average','error','n'};
            else
                % Compute one or more functions
                funs = fieldnames(fun);
                nrFuns = numel(funs);
                fprintf('Applying %d functions (%s) to %d elements\n',nrFuns,strjoin(funs,"/"),size(M,1))
                for f = 1:nrFuns
                    thisFun = funs{f};
                    if isstruct(fun.(thisFun))
                        thisOptions = namedargs2cell(fun.(thisFun));
                    elseif isempty(fun.(thisFun)) || isa(fun.(thisFun),'function_handle')
                        thisOptions = {};
                    else
                        error('Function options for %s must be a struct, empty, or a function handle', thisFun);
                    end
                    data = {M}; % Inputs indexed by group row                  
                    switch thisFun
                        case "average"
                            funN = @(signal,srate) ns.cache.do_average(signal,srate,thisOptions{:});                           
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
                    pool = nsParPool;
                    % Temp cell to store results
                    xCell = cell(nrGrps,1);
                    idv  = repmat("",1, nrGrps);
                    if isempty(pool)
                        for iGrp = 1:nrGrps
                            groupData = cellfun(@(d) selectDataRow(d,iGrp),data,UniformOutput=false);
                            [xCell{iGrp},idv(iGrp)] = funN(groupData{:},srate);
                        end
                    else
                        progressQueue = parallel.pool.DataQueue;                        
                        fprintf('Applying %s in parallel: 0/%d',thisFun,nrGrps);
                        parfor iGrp = 1:nrGrps
                            groupData = cellfun(@(d) selectDataRow(d,iGrp),data,UniformOutput=false);
                            [xCell{iGrp},idv(iGrp)] = funN(groupData{:},srate); %#ok<PFBNS>
                            send(progressQueue,1);
                        end
                        fprintf('\n');
                    end
                    R = vertcat(xCell{:}); % Results table
                    idv = unique(idv); % Should all be the same.                     
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
                        assert(all(size(R{1,col})==size(R{1,idv(1)})),"IDV and DV must match in size")
                    end
                    
                    % Pass the dependent variables to the next computation
                    M = table2cell(R(:,dv));        
                    
                    % Combine the results with G and store idv->dv mapping                    
                    D(idv)  = dv;                     
                    G = [G R]; %#ok<AGROW>                    
                end                       
           end

            % Sort in consistent order - not matched to the tbl query
            G= sortrows(G,intersect(["subject" "session_date" "starttime" "paradigm"  "condition" "channel" "trial"],G.Properties.VariableNames,'stable'));
          
        end
    end



    methods (Static)
        % Compute functions that take a signal with some options and return
        % a table with one or more output columns. Note that each column
        % should contain a row vector of results.
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
