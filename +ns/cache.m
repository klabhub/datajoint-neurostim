classdef (Abstract) cache < handle
    % Abstract superclass used by ns.Epoch and ns.Epoch
    % When working with a set of T/Epochs one often wants to compute
    % derived measures (e.g. a spectrum from a signal), average in
    % different ways (across trials or channels, or subjects), or visualize
    % raw or computed data.
    % This cache class prevents multiple round trips to the server to fetch
    % the data. Instead, the data are fetched once and stored internally in
    % a Matlab table (.T) that also tracks the various primary keys of each
    % row in the DJ table.
    %
    % EXAMPLE
    %
    %
    % BK - Dec 2025
    properties (Constant)
        GROUPVARS = ["subject" "session_date" "starttime" "condition" "trial" "channel"];
    end
    properties (GetAccess =public,SetAccess = protected)
        T (:,:) table  = table;  % The Matlab table that stores the data
        qry (1,1) string =""     % The query that fetched the data
        independent (1,:) string = "time" % The name(s) of the independent variables
        dependent (1,:) string = "signal"  % The name(s) of the dependent variables
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
            G = compute(o,"msten",x=xName,y=yName,average= pv.average);
            x = G{1,xName};
            if xName =="time" && numel(x) ==3
                x = linspace(x(1),x(2),x(3));
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
                    if mod(tileCntr,pv.tilesPerPage)==0
                        if i>1 && pv.linkAxes
                            linkaxes(gcf().Children().Children())
                        end
                        figure;
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
                    imagesc(x,1:nrTrials, G.mean{i})
                    axis xy
                    n = mean(G.n{i},"all");
                    ylabel ("")
                    titlePV= setdiff(["paradigm" rasterGrouping],"",'stable');
                    ttlStr = strjoin(string(G{i,titlePV}),"/");
                else
                    m = G.mean(i,:);
                    ste = G.ste(i,:);
                    n = mean(G.n(i,:));                    
                    h = [h plot(x,m)];                %#ok<AGROW>
                    p = patch([x flip(x)],[m+ste flip(m-ste)],h(end).Color,FaceAlpha= 0.5);
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
                        y = m - matchG.mean(reference,:);
                        ste = ste +matchG.ste(reference,:);
                        h = [h plot(x,y)];                     %#ok<AGROW>
                        p = patch([x flip(x)],[y+ste flip(y-ste)],h(end).Color,FaceAlpha= 0.5);
                        p.EdgeColor = h(end).Color;
                        legStr = [legStr G.condition(i)+"-"+ pv.delta]; %#ok<AGROW>
                    end
                end
            end
            if ~isempty(legStr); legend(h,legStr); end % For the last tile
            if pv.linkAxes
                linkaxes(gcf().Children().Children())
            end
        end

        function [G,dv,idv] =  compute(o,fun,options,pv)
            % Compute derived measures from the EpochChannel table.
            % fun - fft, psd, pmtm, snr, msten
            %       fft - Uses fft() to compute amplitude, phase as a function of frequency
            %        psd - Uses pspectrum to compute power spectral density as a function of frequency
            %        pmtm - Uses multittaper pmtm to compute a spectrogram as a function of time and frequency
            % options - Struct passed to the compute function - different
            %           for each fun.
            %        fft - no options
            %        psd
            %       pmtm
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
                fun  (1,1) string {mustBeMember(fun,["fft" "psd" "pmtm" "snr" "msten" "wavelet"])}
                options (1,:) struct = struct([]);  % Options for the fun
                pv.channel (:,1) double = []            % Select a subset of channels
                pv.trial (:,1) double = []            % Select a subset of trials
                pv.timeWindow (1,2) double = [-inf inf]  % Select a time window to operate on
                pv.average (1,:) string {mustBeMember(pv.average,["" "subject" "session_date" "starttime" "condition" "trial" "channel"])} = ["trial" "channel"]
                pv.x (1,1) string = o.independent  % Which  column in .T to use on the horizontal axis
                pv.y (1,1) string = o.dependent    % Column in .T to use as the dependent variable
            end
            fill(o);% Fill the cache

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
            if any(isfinite(pv.timeWindow))
                assert(ismember("time",o.T.Properties.VariableNames),"timeWindow restriction can only be used on a cache with a time column.")
                % Crop to the timeWindow for this operation.
                t = restrictedT{1,"time"}; % Time in secods (Taking first row as all should be the same)
                if numel(t)==3
                    t = linspace(t(1),t(2),t(3));
                end
                keep = do.ifwithin(t,pv.timeWindow/1000);
                restrictedT(:,pv.y)= rowfun(@(x) {x(keep,:)},restrictedT(:,pv.y) ,'ExtractCellContents',true);
                t= t(keep);
                assert(~isempty(t),'No time points left in the analysis window ([%f %f])',pv.timeWindow(1),pv.timeWindow(2));
                restrictedT.time = repmat([t(1) t(end) numel(t)],height(restrictedT),1);
            end
            if isempty(restrictedT)
                error('No data in this table');
            end

            %%  Average/group
            if pv.average~=""
                grouping = setdiff(ns.cache.GROUPVARS,pv.average,'stable');
                [grp,G] = findgroups(restrictedT(:,grouping));
                % Average per group
                if fun=="msten"
                    M = splitapply(@ns.cache.msten,restrictedT.(pv.y),grp);
                    
                else
                    if iscell(restrictedT{1,pv.y})
                        M = splitapply(@(x) {mean(cat(2,x{:}),2,"omitmissing")},restrictedT.(pv.y),grp);
                    else
                        M = splitapply(@(x) {mean(x,1,"omitmissing")'},restrictedT.(pv.y),grp);
                    end
                end
            else
                G = restrictedT;
                M = restrictedT;
            end
            nrGrps = height(M);

            %% Determine which function to compute
            % Map string to function handle and do error checking
            if fun=="msten"
                dv = "mean";
                idv = "time";
                % Combine with G
                G = [G M];
            else
                switch fun
                    case "fft"
                        assert(isempty(options),"fft does not take any options")
                        fun = @(x) ns.cache.do_fft(x,o.samplingRate);
                        dv = ["amplitude" "phase"];
                        idv = "frequency";
                    case "psd"
                        if isempty(options)
                            opts ={};
                        else
                            opts = namedargs2cell(options);
                        end
                        fun = @(x) ns.cache.do_psd(x,o.samplingRate,opts{:});
                        idv = "frequency";
                        dv = "power";
                    case "snr"
                        fun = @(x) ns.cache.do_snr(x,o.samplingRate,options);
                        idv = "frequency";
                        dv = "power";
                    case "pmtm"
                        if ~isfield(options,"args")
                            % Defaults
                            args = {4,7}; % time-halfbandwidth product =4, averaging weights =7
                        else
                            args = options.args;
                        end
                        options = rmfield(options,'args');
                        if isempty(options)
                            options = {};
                        else
                            options = namedargs2cell(options); % The rest as parm/value pairs
                        end
                        fun = @(x) ns.cache.do_pmtm(x,args{:},o.samplingRate,options{:});
                        idv = "frequency";
                        dv = "power";
                    case "wavelet"
                        if isempty(options)
                            options= struct('limits', [0.5 100],'fwhm',[2 0.2],'nfrex',40);
                        end
                        fun = @(x) ns.cache.do_wavelet(x,o.samplingRate,options.fwhm,options.nfrex,options.limits);
                        idv = ["frequency" "time"];
                        dv = "power";
                    otherwise
                        error('Unknown function %s', fun);
                end

                %% Apply the fun to the mean signal
                results = splitapply(fun,M,(1:nrGrps)');
                % Combine with G
                G = [G results];               
            end
            % Combine with align/time/paradigm information. Note this
                % assumes these are constant across the group (picking
                % only the first here). fill() assures this is the case.
                if pv.average ~=""
                    P = groupsummary(restrictedT, grouping, @(x) x(1,:), ["align" pv.x "paradigm"]);
                    P = renamevars(P,["fun1_align" "fun1_"+pv.x "fun1_paradigm"],["align" pv.x "paradigm"]);
                    G =innerjoin(G,P);
                end
            % Sort in consistent order - not matched to the tbl query
            G= sortrows(G,intersect(["subject" "session_date" "starttime" "paradigm"  "condition" "channel" "trial"],G.Properties.VariableNames,'stable'));
        end
    end



    methods (Static)
        % Compute functions that take a signal with some options and return
        % a table with one or more output columns. Note that each column
        % should contain a row vector of results.
        function v = do_fft(signal, fs)
            % do_fft - Computes FFT amplitude and phase for each
            %               epoch. Only includes real frequencies.
            %
            % Outputs (table columns):
            %   amplitude: Amplitude of the FFT.
            %   phase: Phase of the FFT.
            %   frequency: Corresponding real frequencies.

            signal =cat(2,signal{:}); % Concatenate epochs
            % Compute FFT for each slice along time dim 1
            fftResult = fft(signal);

            % Calculate real frequencies
            N = size(signal, 1);
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
            v= table(amplitude(:)',phase(:)',freq(:)','VariableNames',{'amplitude','phase','frequency'});
        end
        function v = do_psd(signal, fs, varargin)
            % Power spectral density.
            % Table with power and frequency
            signal =cat(2,signal{:}); % Concatenate epochs
            signal = signal - mean(signal,1,"omitmissing");
            [power, freq] = pspectrum(signal, fs, varargin{:});
            v= table(power(:)',freq(:)','VariableNames',{'power','frequency'});
        end
        function v = do_pmtm(signal, varargin)
            % Multitaper power and frequency
            signal =cat(2,signal{:}); % Concatenate epochs
            [power, freq] = pmtm(signal, varargin{:});
            % Make table, force rows
            v = table(power(:)',freq(:)','VariableNames',{'power','frequency'});
        end
        function v = do_wavelet(signal,fs, fwhm,nfrex,limits)
            % Code adapted from Cohen M. X. (2019). A better way to
            % define and describe Morlet wavelets for time-frequency
            % analysis. NeuroImage, 199, 81-86.
            % https://doi.org/10.1016/j.neuroimage.2019.05.048
            signal =cat(2,signal{:}); % Concatenate epochs
            nrSamples= size(signal,1);
            % time-frequency parameters
            freq  = linspace(limits(1),limits(2),nfrex)';
            fwhm = linspace(fwhm(1),fwhm(2),nfrex)'; % variable fwhm
            assert(all(fwhm.*freq>=1),"The FWHM is too small (should have more than one cycle per window)");

            % setup wavelet and convolution parameters
            wavet = (-5:1/fs:5)';
            halfw = floor(length(wavet)/2)+1;
            nConv = nrSamples + length(wavet) - 1;

            % initialize time-frequency matrix
            spectrogram = zeros(nfrex,nrSamples);

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
            v = table({power},freq(:)','VariableNames',{'power','frequency'});
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

        function v = msten(x)
            X =cat(2,x{:});
            v = {mean(X,2,"omitmissing")', ...  % Mean
                (std(X,0,2,"omitmissing")./sqrt(sum(~isnan(X),2,"omitmissing")))',...  % Standard error
                sum(~isnan(X),2,"omitmissing")'};  % Non-Nan N
            % Make a table.
            v = cell2table(v,"VariableNames",{'mean','ste','n'});
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
                assert(isscalar(unique(preFetch.time(:,3))),'Rows of the EpochChannel table must have the same numbers of samples.');
                assert(isscalar(unique({preFetch.align.plugin})),'Rows of the EpochChannel should be aligned to the same plugin.');
                assert(isscalar(unique({preFetch.align.event})),'Rows of the EpochChannel should be aligned to the same event.');
                o.T =fetchtable(src,'*','ORDER BY channel');
                o.qry = src.sql;
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