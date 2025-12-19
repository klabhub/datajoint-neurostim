classdef (Abstract) cache < handle

    properties (GetAccess =public,SetAccess = protected)  
        T (:,:) table  = table;
        qry (1,1) string =""  % Stores the query that fetched the data     
    end
   

    methods (Access = public)
      function plot(o,pv)
            % Plot evoked responses for all rows in the table.
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
                pv.raster (1,1) logical = false         % Set to true to show trials as rasters (removes "trial" from pv.average)
            end
            fill(o);% Fill the cache

            % Start a new tile when these values change:
            newTileEach = ["paradigm" "subject" "session_date" "starttime"];
            dimension = unique(o.T.dimension);
            %% Fill the cache, then perform averaging per group
            
            if pv.raster
                pv.average = setdiff(pv.average,"trial");
            end
            % Determine mean, ste and n averaging over the pv.average
            % dimension (anything that is left out of grouping)
            grouping = setdiff(["subject" "session_date" "starttime" "paradigm"  "condition" "trial" "channel"  ],pv.average,'stable');            
            G = ns.cache.cacheCompute(o.T,@ns.EpochChannel.do_msten,["mean" "ste" "n"],channel=pv.channel,trial=pv.trial,grouping=grouping);
            
            if pv.raster
                % Concatenate the trials into a raster matrix in G.
                grouping = setdiff(grouping,"trial",'stable');
                P = groupsummary(G, grouping, @(x) x(1,:), ["align" "time"]);
                P = renamevars(P,["fun1_align" "fun1_time" ],["align" "time"]);
                G  = groupsummary(G,grouping,@(x) ({cat(1,x)}),["mean" "ste" "n"]);                
                G= renamevars(G,["fun1_mean" "fun1_ste" "fun1_n"],["mean" "ste" "n"]);                                       
                G =innerjoin(G,P);
                newTileEach = union(newTileEach,"condition");          
                G= sortrows(G,intersect([ "subject" "session_date" "starttime" "condition" "channel" "trial" "paradigm"],G.Properties.VariableNames,'stable'));
            end
                   
            
            %% Figure
            tileCntr=0;
            nrTimeSeries = height(G);
            newTileEach = intersect(newTileEach,G.Properties.VariableNames);
            legStr =string([]);
            for i = 1:nrTimeSeries
                t =  G.time(i,:);
                t = linspace(t(1),t(2),t(3));
                align =G.align(i);                
                if i==1 || (~isempty(newTileEach) && any(G{i,newTileEach} ~= G{i-1,newTileEach}))
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
                if pv.raster
                    % Show each condition in a separate tile
                    nrTrials= size(G.mean{i},1);
                    imagesc(t,1:nrTrials, G.mean{i})
                    n = mean(G.n{i},"all");
                    ylabel ("Trial")
                    titlePV= setdiff(["paradigm" grouping],"",'stable');
                    ttlStr = strjoin(string(G{i,titlePV}),"/");
                else
                    m = G.mean(i,:);
                    ste = G.ste(i,:);
                    n = mean(G.n(i,:));
                    h = [h plot(t,m)];                %#ok<AGROW>
                    p = patch([t flip(t)],[m+ste flip(m-ste)],h(end).Color,FaceAlpha= 0.5);                
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
                        h = [h plot(t,y)];                     %#ok<AGROW>
                        p = patch([t flip(t)],[y+ste flip(y-ste)],h(end).Color,FaceAlpha= 0.5);
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


        function [T,dv,iv] = compute(o,fun, options,pv)
            % Compute derived measures from the EpochChannel table.
            % fun - fft, psd, pmtm, snr, msten
            % options - Struct passed to the compute function - different
            %           for each fun
            % channel  - Select a subset of channels
            % trial    - Select a subset of trials
            % grouping - Group analysis by these fields. Signals are
            %            averaged within group before applying the fun.
            arguments
                o (1,1) 
                fun  (1,1) string {mustBeMember(fun,["fft" "psd" "pmtm" "snr" "msten"])}
                options (1,:) struct = struct([]);  % Options for the fun
                pv.channel (:,1) double = []            % Select a subset of channels
                pv.trial (:,1) double = []            % Select a subset of trials
                pv.timeWindow (1,2) double = [-inf inf]  % Select a time window to operate on
                pv.grouping (1,:) string {mustBeMember(pv.grouping,["subject" "session_date" "starttime" "condition" "trial" "channel"])} = ["subject" "session_date" "starttime" "condition"]
            end           
            fill(o);% Fill the cache            
            samplingRate = o.samplingRate;
            switch fun
                case "msten"
                    fun =@ns.cache.msten;                    
                    names = ["mean" "ste" "n"];
                    dv   = names;
                    iv   = "";  
                case "fft"
                    assert(isempty(options),"fft does not take any options")
                    fun = @(x) ns.cache.do_fft(x,samplingRate);
                    names = ["fft" "phase" "frequency"];
                    dv   = "fft";
                    iv   = "frequency"; 
                case "psd"
                   if isempty(options)
                        opts ={};
                   else
                        opts = namedargs2cell(options);
                   end
                    fun = @(x) ns.cache.do_psd(x,samplingRate,opts{:});
                    names = ["power" "frequency"];
                    dv   = "power";
                    iv   = "frequency"; 
                case "pmtm"                   
                    % opts.args = {4,7}
                    % opts.Tapers= 'Sine'
                    assert(isfield(options,"args"),"pmtm needs an .args field in the options struct, which will be passsed verbatim to pmtm.")
                    args =options.args;
                    opts = namedargs2cell(rmfield(options,'args')); % The rest as parm/value pairs
                    fun = @(x) ns.cache.do_pmtm(x,args{:},opts{:},samplingRate);
                    names = ["power" "frequency"];
                    dv   = "power";
                    iv   = "frequency"; 
                case "snr"
                    % First compute power using FFT
                    signalFun = @(x) ns.cache.do_fft(x,samplingRate);
                    signalNames = ["fft" "phase" "frequency"];
                    pwrT = ns.cache.cacheCompute(o.T,signalFun,signalNames ,channel =pv.channel,trial =pv.trial,timeWindow = pv.timeWindow,grouping=pv.grouping);

                    fq_res = mode(diff(pwrT{1,"frequency"}));
                    n_bins = round(options.bin/fq_res);
                    n_bins_skip = round(options.bin_skip/fq_res);
                    mid_bin = n_bins + 1;
                    kernel = true(1, 2*n_bins + 1);
                    kernel(mid_bin-n_bins_skip : mid_bin+n_bins_skip) = false;
                    noiseFun = @(x) ns.cache.do_noise(x,kernel);
                    noiseNames = "noise";
                    N = ns.cache.cacheCompute(pwrT,noiseFun,noiseNames,signal="fft",grouping=pv.grouping);

                    snr = cellfun(@(sig,nois) abs(sig)./nois,pwrT.fft,N.noise,'UniformOutput',false);
                    pwrT= addvars(pwrT,snr);
                    pwrT= removevars(pwrT,signalNames);
                    T = pwrT;
                    dv   = "snr";
                    iv   = "frequency"; 
                    return; % All done.
                otherwise 
                    error('Unknown function %s', fun);
            end
            % Send to cacheCompute
            T = ns.cache.cacheCompute(o.T,fun,names,channel =pv.channel,trial =pv.trial,timeWindow = pv.timeWindow,grouping=pv.grouping);
        end
    end

   

    methods (Static)
     
       
        function G = cacheCompute(o,fun, names,pv)
            % This function evaluates a function on cached epochChannel
            % table. It is called internally by plot and by compute.
            arguments
                o table
                fun (1,1)
                names (1,:) string
                pv.channel (:,1) double = []  % Select a subset of channels
                pv.trial   (:,1) double = []  % Select a subset of trials
                pv.timeWindow (1,2) double = [-inf inf] % Select a time window
                pv.grouping (1,: ) string = ["subject" "session_date" "starttime" "condition"]
                pv.signal (1,1) string = "signal"; % Which column of the cached table to operate on
            end
            nrNames  = numel(names);
            if ~isempty(pv.channel)
                % Restrict channels for this operation
                o = o(ismember(o.channel,pv.channel),:);
            end
            if ~isempty(pv.trial)
                % Restrict trials for this operation
                o = o(ismember(o.trial,pv.trial),:);
            end
            if any(isfinite(pv.timeWindow))
                % Crop to the timeWindow for this operation.
                t = o{1,"time"}; % Time in secods
                t = linspace(t(1),t(2),t(3));
                keep = do.ifwithin(t,pv.timeWindow/1000);
                o(:,pv.signal)= rowfun(@(x) {x(keep,:)},o(:,pv.signal) ,'ExtractCellContents',true);
                t= t(keep);
                assert(~isempty(t),'No time points left in the analysis window ([%f %f])',pv.timeWindow(1),pv.timeWindow(2));
                o.time = repmat([t(1) t(end) numel(t)],height(o),1);
            end
            if isempty(o)
                error('No data in this table');
            end
            
            [grp,G] = findgroups(o(:,pv.grouping));
            % Average first within group
            M = splitapply(@(x) {mean(cat(2,x{:}),2)},o.(pv.signal),grp);
            nrGrps = height(M);
            % Then apply the fun to the mean signal
            results = splitapply(fun,M,(1:nrGrps)');

            assert(size(results,2)==nrNames,"The fun (%s) returns %d values, but %d names have been provided",func2str(fun),size(results,2),nrNames);
            % Rename
            for r=1:nrNames
                sz =cellfun(@size,results(:,r),'uni',false);
                sz = cat(1,sz{:});
                sz = unique(sz,"rows");
                canBeDouble = size(sz,1)==1 && sz(1)==1;
                if canBeDouble
                    G= addvars(G, cat(1,results{:,r}),'NewVariableNames',names(r));
                else
                    G = addvars(G, results(:,r), 'NewVariableNames', names(r));
                end
            end
            % Combine with align/time/paradigm information. Note this
            % assumes these are constant across the group (picking
            % only the first here). fillCache assures this is the case.
            
                
                P = groupsummary(o, pv.grouping, @(x) x(1,:), ["align" "time"]);
                P = renamevars(P,["fun1_align" "fun1_time" ],["align" "time"]);
                G =innerjoin(G,P);
            
            % Sort in consistent order - not matched to the tbl query
            G= sortrows(G,intersect(["subject" "session_date" "starttime" "paradigm"  "condition" "channel" "trial"],G.Properties.VariableNames,'stable'));
        end

        function v = do_fft(signal, fs)
            % fftSliceReal - Computes FFT amplitude and phase for each
            %               epoch. Only includes real frequencies.
            %
            % Inputs:
            %   data: Input multi-dimensional data.
            %   fs: Sampling frequency
            %
            % Outputs:
            %   amplitude: Amplitude of the FFT.
            %   phase: Phase of the FFT.
            %   frequencies: Corresponding real frequencies.

            signal =cat(2,signal{:}); % Concatenate epochs
            % Compute FFT for each slice along time dim 1
            fftResult = fft(signal);

            % Calculate real frequencies
            N = size(signal, 1);
            if mod(N, 2) == 0
                frequencies = (0:N/2) * fs / N;
                idx = 1:N/2+1;
            else
                frequencies = (0:(N-1)/2) * fs / N;
                idx = 1:(N+1)/2;
            end

            amplitude = 2*abs(fftResult(idx,:,:)/sqrt(N));
            phase = angle(fftResult(idx,:,:));
            % Return as single cell for cacheCompute
            v= {amplitude,phase,frequencies};
        end
        function v = do_psd(signal, fs, varargin)
            % Power spectral density
            signal =cat(2,signal{:}); % Concatenate epochs
            [power, frequencies] = pspectrum(signal, fs, varargin{:});
            v= {power,frequencies};
        end
        function v = do_pmtm(signal, varargin)
            signal =cat(2,signal{:}); % Concatenate epochs
            [power, frequencies] = pmtm(signal, varargin{:});
            v = {power,frequencies};
        end
        function v = do_wavelet_sum(signal, fs, freq_range,wid_range, n_frex, n_wid)

            % Code adapted from Nicole, then from Cohen M. X. (2019). A better way to
            % define and describe Morlet wavelets for time-frequency
            % analysis. NeuroImage, 199, 81-86.
            % https://doi.org/10.1016/j.neuroimage.2019.05.048

            arguments

                signal (:,:,:)
                fs (1,1) double
                freq_range (1,2) {mustBeNumeric, mustBePositive, mustBeInteger}
                wid_range (1,2) {mustBeNumeric, mustBePositive, mustBeInteger}
                n_frex (1,1) {mustBeNumeric, mustBePositive, mustBeInteger} = 80
                n_wid (1,1) {mustBeNumeric, mustBePositive, mustBeInteger} = 20



            end
            %%% time-frequency parameters
            frex  = linspace(freq_range(1), freq_range(2),n_frex);
            fwhm = linspace(wid_range(1), wid_range(2),n_wid); % variable fwhm


            % setup wavelet and convolution parameters
            wavet = -5:1/fs:5;
            halfw = floor(length(wavet)/2)+1;
            nConv = length(signal) + length(wavet) - 1;

            % initialize time-frequency matrix
            tf = zeros(length(frex),length(tidx));


            % spectrum of data
            dataX = fft(signal,nConv);

            % loop over frequencies
            for fi=1:length(frex)

                % create wavelet
                waveX = fft( exp(2*1i*pi*frex(fi)*wavet).*exp(-4*log(2)*wavet.^2/fwhm(fi).^2),nConv );
                waveX = waveX./max(waveX); % normalize
                % convolve
                as = ifft( waveX.*dataX );
                % trim
                as = as(halfw:end-halfw+1);

                % power
                power = abs(as).^2;

                % storing power results in time-frequency matrix;
                % frequency x time
                tf(fi,:) = 10*log10(power(tidx)); % no baseline normalization - used for resting-state

                % empirical FWHM
                hz = linspace(0,freqSamp,nConv);
                idx = dsearchn(hz',frex(fi));
                fx = abs(waveX);
                empfwhm(4,fi) = hz(idx-1+dsearchn(fx(idx:end)',.5)) - hz(dsearchn(fx(1:idx)',.5));

                % s  = fwhm*(2*pi-1)/(4*pi); % normalized width
                % x  = hz-peakf;             % shifted frequencies
                % fx = exp(-.5*(x/s).^2);    % gaussian
            end
            % average power
            v = mean(tf,2);
        end
        function v= do_noise(signal,kernel)
            % Mean power in neighboring frequency bins.
            kernel = kernel(:); % Force column to operate along samples
            signal =cat(2,signal{:}); % Concatenate epochs
            signal = abs(signal);
            kernelWidth = floor(numel(kernel)/2);
            noise = conv2(signal, kernel, 'same');% Uses zero padding
            noise(1:kernelWidth,:) = NaN;
            noise(end-kernelWidth+1:end,:) = NaN;
            noise =noise/sum(kernel); % Average instead of sum
            v = {noise};
        end
        function v = do_msten(x)
            % Determine mean, standard error, and N and return as a single cell.
            % Note that time dimension must be rows for cacheCompute to
            % work
            X =cat(2,x{:});
            v = {mean(X,2,"omitmissing")', ...  % Mean
                (std(X,0,2,"omitmissing")./sqrt(sum(~isnan(X),2,"omitmissing")))',...  % Standard error
                sum(~isnan(X),2,"omitmissing")'};  % Non-Nan N
        end
    end



    methods (Access= protected)
        function fill(o)
            % Fetch the data if the underlying query has changed                         
            src = getCacheQuery(o);
            if string(src.sql) ~= o.qry
                newT=fetchtable(src,'*','ORDER BY channel');
                % Safety check; time and align should match for all rows
                % in the table.
                t= newT.time;
                assert(isscalar(unique(t(:,3))),'Rows of the (Te/E)pochChannel table must have the same numbers of samples.')
                assert(isscalar(unique({newT.align.plugin})),'Rows of the (Te/E)pochChannel should be aligned to the same plugin')
                assert(isscalar(unique({newT.align.event})),'Rows of the (Te/E)pochChannel should be aligned to the same event')
                o.T = newT;
                o.qry = src.sql;
            end
        end
    end
    
    methods (Abstract, Access = protected)
        src = getCacheQuery(o)
    end
end