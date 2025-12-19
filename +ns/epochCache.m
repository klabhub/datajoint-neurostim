classdef epochCache   <handle

    properties (GetAccess =public,SetAccess = protected)  
        T (:,:) table  = table;
        cacheQry (1,1) string =""  % Stores the query that fetched the data     
    end

    

    methods 
     
       
        function G = compute(o,fun, names,pv)
            % This function evaluates a function on cached epochChannel
            % table. It is called internally by plot and by compute.
            arguments
                o epochCache
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
                o = o(ismember(o.T.channel,pv.channel),:);
            end
            if ~isempty(pv.trial)
                % Restrict trials for this operation
                o = o(ismember(o.T.trial,pv.trial),:);
            end
            if any(isfinite(pv.timeWindow))
                % Crop to the timeWindow for this operation.
                t = o.T{1,"time"};
                t = linspace(t(1),t(2),t(3));
                keep = do.ifwithin(t,pv.timeWindow);
                o(:,pv.signal)= rowfun(@(x) {x(keep,:)},o(:,pv.signal) ,'ExtractCellContents',true);
                t= t(keep);
                o.T.time = repmat([t(1) t(end) numel(t)],height(o),1);
            end
            if isempty(o.T)
                error('No data in this table');
            end
            [grp,G] = findgroups(o.T(:,pv.grouping));
            % Average first within group
            M = splitapply(@(x) {mean(cat(2,x{:}),2)},o.T.(pv.signal),grp);
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
            P = groupsummary(o.T, pv.grouping, @(x) x(1,:), ["align" "time"]);
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



    methods (Access={?ns.EpochChannel,?nsTepochChannel} )
        function fill(o,relvar)
            % Fetch the data if the underlying query has changed                         
            if string(relvar.sql) ~= o.cacheQry
                newT=fetchtable(relvar,'paradigm','signal','time','condition','align','ORDER BY channel');
                % Safety check; time and align should match for all rows
                % in the table.
                t= newT.time;
                assert(isscalar(unique(t(:,3))),'Rows of the (Te/E)pochChannel table must have the same numbers of samples.')
                assert(isscalar(unique({newT.align.plugin})),'Rows of the (Te/E)pochChannel should be aligned to the same plugin')
                assert(isscalar(unique({newT.align.event})),'Rows of the (Te/E)pochChannel should be aligned to the same event')
                o.T = newT;
                o.cacheQry = relvar.sql;
            end
        end
    end

end