%{
# Preprocessed and epoched data
-> ns.Epoch
channel : int
---
signal : longblob
%}

classdef EpochChannel < dj.Part & dj.DJInstance

    properties
        data
        
        timepoints_ = []
        frequencies_ = []

    end

    properties (SetAccess = protected)

        master = ns.Epoch

    end

    properties (Dependent)

        Epoch
        C
        timepoints % time points associated with ecTbl.data
        channels
        frequencies
                
    end

    methods
         function varargout = transform(ecTbl, varargin)
                                
                cTbl = ecTbl.C;
                epTbl = ns.EpochParm & ecTbl;
                etags = unique(epTbl{'etag'},'rows');
                ctags = unique(cTbl{'ctag'},'rows');
                files = unique(cTbl{'filename'}, 'rows');


                if size(etags,1)>1 || size(ctags,1)>1 || size(files,1) > 1
                    error("Transformed data contains data from different files and/or preprocessing pipelines.");
                end

                if isempty(ecTbl.data)

                    ep = ecTbl.load();

                else

                    ep = ecTbl.data;

                end
                n_arg = nargin - 1;
                ii = 1;

                while ii <= n_arg

                    argN = varargin{ii};
                    ii = ii + 1;

                    switch argN
                        case 'include'

                            argN = varargin{ii};
                            ep = ecTbl.include(ep, argN{:});
                            
                            ii = ii + 1;
                            if ii > n_arg

                                varargout{1} = ep;

                            end
                            
                        case 'exclude'
                        case 'baseline'
                        case 'detrend'
                        case 'fft'

                            [ep, ph, fq] = ecTbl.do_fft(ep, cTbl.srate, 3);
                            ecTbl.frequencies_ = fq;
                            if ii > n_arg

                                varargout{1} = ep;
                                varargout{2} = ph;
                                varargout{3} = fq;

                            end
                        case 'psd'

                            if ii <= n_arg
                                opts = varargin{ii};
                                if iscell(opts)
                                    ii = ii + 1;
                                else
                                    opts = {};
                                end
                            end
                            [ep, fq] = ecTbl.do_psd(ep, cTbl.srate, opts{:});

                            if ii > n_arg

                                varargout{1} = ep;
                                varargout{2} = fq;

                            end

                            ecTbl.frequencies_ = fq;

                        case 'pmtm'

                            if ii <= n_arg
                                opts = varargin{ii};
                                if iscell(opts)
                                    ii = ii + 1;
                                else
                                    opts = {};
                                end
                            end
                            [ep, fq] = ecTbl.do_pmtm(ep, opts{:});

                            if ii > n_arg

                                varargout{1} = ep;
                                varargout{2} = fq;

                            end

                            ecTbl.frequencies_ = fq;


                        case 'snr'    

                            argN = varargin{ii}; %options
                            ii = ii + 1;

                            fq_res = ecTbl.frequencies(2)-ecTbl.frequencies(1);
                            fq_bin_skip = argN{1};
                            fq_bin = argN{2};
                            
                            n_bins = round(fq_bin/fq_res);
                            n_bins_skip = round(fq_bin_skip/fq_res);
                            mid_bin = n_bins + 1;
                            kernel = true(1, 2*n_bins + 1);
                            kernel(mid_bin-n_bins_skip : mid_bin+n_bins_skip) = false;

                            if numel(argN) > 2
                                method = argN{3};
                            else
                                method = @mean;
                            end

                            noise = ecTbl.calculate_noise_(ep, kernel, 3, method);
                            ep = ep./noise;

                            if ii > n_arg
                                varargout{1} = ep;
                                varargout{2} = ecTbl.frequencies;
                            end

                    end

                end

                ecTbl.data = ep;
                
        end

        function [ep, trl, ch] = load(ecTbl, qry)

            % loads data
            if nargin > 1 && ~isempty(qry)

                ecTbl = ecTbl & qry;          

            end
            %% Fetch data                
            t_fetch = tic;
            fprintf('Fetching epochs...');
            % Change it so that ep, trl, and ch is separate for every ctag
            % and etag, currently it should only load when that is the case
            dat = fetch(ecTbl , 'signal');
            ep = cat(1, dat(:).signal);
            trl = cat(1, dat(:).trial);
            ch = cat(1, dat(:).channel);

            % get dimension size of the final ep matrix
            n_ep = numel(unique(trl));
            n_ch = numel(unique(ch));
            n_t = size(ep,2);

            % reshape ep to n_ch by n_ep by n_t
            ep = reshape(ep, n_ch, n_ep, n_t);
            %permute to get the desired shape n_ep n_ch, n_t
            ep = permute(ep, [2, 1, 3]);

            ecTbl.data = ep;

            fprintf("completed.\n"); toc(t_fetch);

        end

        function ep = include(ecTbl, ep, pv)

            arguments

                ecTbl ns.EpochChannel
                ep double
                pv.time_window = []
                pv.frequency_window = []
                
                pv.trials = []
                pv.conditions = []
            
            end            

            eTbl = ecTbl.Epoch;

            if ~isempty(pv.time_window)

                isTIn = do.ifwithin(ecTbl.timepoints, pv.time_window);
                ecTbl.timepoints_ = ecTbl.timepoints(isTIn);
                ep = ep(:,:,isTIn);

            end

            if ~isempty(pv.frequency_window)

                isTIn = do.ifwithin(ecTbl.frequencies, pv.frequencies_window);
                ecTbl.frequencies_ = ecTbl.frequencies(isTIn);
                ep = ep(:,:,isTIn);

            end

            if ~isempty(pv.trials)

                trl = fetch(eTbl,'trial');
                trl = [trl(:).trial];
                isTrlIn = ismember(trl, pv.trials);

                ep = ep(isTrlIn, :, :);

            end

            if ~isempty(pv.conditions)

                isTrlIn = ismember(eTbl.conditions, pv.conditions);
                ep = ep(isTrlIn, :, :);

            end      
        end % include()

        
    end

    methods %Get methods

        function eTbl = get.Epoch(self)

            eTbl = ns.Epoch & self;

        end

        function cTbl = get.C(self)

            cTbl = ns.C & self;

        end

        function t = get.timepoints(self)

            if isempty(self.timepoints_)

                t = self.Epoch.timepoints;

            else

                t = self.timepoints_;

            end

        end

        function f = get.frequencies(self)

            f = self.frequencies_;

        end

        function ch = get.channels(self)

            ch = fetch(self, 'channel');
            ch = unique([ch(:).channel]');

        end

        
    end

    methods (Static)

        function [amplitude, phase, frequencies] = do_fft(data, fs, n)
            % fftSliceReal - Computes FFT amplitude and phase for each slice along the n'th dimension.
            %                Only includes real frequencies.
            %
            % Inputs:
            %   data: Input multi-dimensional data.
            %   n: Dimension along which to slice and compute FFT.
            %   fs: Sampling frequency (optional, default is 1).
            %
            % Outputs:
            %   amplitude: Amplitude of the FFT.
            %   phase: Phase of the FFT.
            %   frequencies: Corresponding real frequencies.

            n_dim = ndims(data);

            % Compute FFT for each slice
            fftResult = fft(data, [], n_dim);

            idx = repmat({':'},1,n_dim);
            % Calculate real frequencies
            N = size(data, n);
            if mod(N, 2) == 0
                frequencies = (0:N/2) * fs / N;
                idx{n_dim} = 1:N/2+1;
            else
                frequencies = (0:(N-1)/2) * fs / N;
                idx{n_dim} = 1:(N+1)/2;
            end

            amplitude = 2*abs(fftResult(idx{:})/sqrt(N));
            phase = angle(fftResult(idx{:}));


        end

        function [p, fq] = do_psd(data, fs, varargin)

            n_ep = size(data,1);
            
            p = cell(n_ep,1);
            for ii = 1:n_ep

                dataN = squeeze(data(ii,:,:))';
                [pN, fq] = pspectrum(dataN, fs, varargin{:});
                p{ii} = pN';

            end

            p = permute(cat(ndims(data),p{:}),[3,1,2]);
        end

        function [p, fq] = do_pmtm(data, varargin)

            n_ep = size(data,1);
            
            p = cell(n_ep,1);
            for ii = 1:n_ep

                dataN = squeeze(data(ii,:,:))';
                [pN, fq] = pmtm(dataN, varargin{:});
                p{ii} = pN';

            end

            p = permute(cat(ndims(data),p{:}),[3,1,2]);

        end

        function [] = do_wavelet_sum(signal, fs, freq_range,wid_range, n_frex, n_wid)

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
            % initialize power matrix
            avgPower = zeros(length(frex));

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
            mean_power = mean(tf,2);

        end

        function noise_arr = calculate_noise_(arr, filter_kernel, n_dim, method)
            
            arr = max(arr, eps);
            if nargin < 3 || isempty(n_dim), n_dim = ndims(arr); end
            if nargin < 4; method = @mean; end

            idx_selector = repmat({':'}, 1, ndims(arr));
            idx_assigner = repmat({':'}, 1, ndims(arr));
            
            n_points = size(arr, n_dim);
            noise_arr = zeros(size(arr));
            n_kernel = numel(filter_kernel);

            for i = 1:n_points
                idx_assigner{n_dim} = i;

                if i <= n_kernel/2 || i > n_points - n_kernel/2

                    noise_arr(idx_assigner{:}) = NaN;

                else

                    idxN = false(1, size(arr,n_dim));
                    idxN(floor(i+1-n_kernel/2):floor(i+n_kernel/2)) = filter_kernel;

                    idx_selector{n_dim} = idxN;

                    noise_arr(idx_assigner{:}) = method(arr(idx_selector{:}), n_dim);

                end

            end
        end
    end
    
end