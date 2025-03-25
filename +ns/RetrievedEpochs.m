classdef RetrievedEpochs < handle

    %{
A subset of epochs retrieved / fetched from ns.Epochs for plotting,
transforming data etc.
The retrieved epochs must call for the same timepoints and channels.
    %}

    properties

        Epoch (1,1) ns.Epoch
        EpochParameter (1,1) ns.EpochParameter
        data table
        channels (1,:) double = []
        time_win double = []        
        trials = {}

        frequencies

    end

    properties (Dependent)

        timepoints
        epoch_win
        dt

    end

    properties (Hidden)

        isPadded_ = false
        pad_filler_ = 'zeros';
        signal_length_after_padding_
        transform_steps_ = {}

    end

    properties (Hidden, Constant)

        separate_figure_by_levels_ = ["ctag", "etag", "dimension", "subject"]

    end

    methods

        function ep = RetrievedEpochs(eTbl, pv)

            arguments

                eTbl (1,1) ns.Epoch

                pv.channels = [] % if empty all channels
                pv.trials = {} % if empty, all trials
                pv.time_window = [];

            end

            ep.Epoch = eTbl;
            ep.EpochParameter = ns.EpochParameter & eTbl;
            if isempty(pv.time_window)

                pv.time_window = ep.EpochParameter.epoch_win;

            end
            
            assert(ep.EpochParameter.nrows == 1, ...
                "retrieve(ns.Epoch,...) is designed for compatible epoch structures corresponding to a single etag.")

            ep.channels = pv.channels;
            ep.time_win = pv.time_window;
            ep.trials = pv.trials;
            ep.retrieve;


        end

        function op = transform(ep, varargin)

            
            op = ep.data;

            iOper = 1;
            while iOper <= nargin-1

                argN = varargin{iOper};
                iOper = iOper + 1;
                switch lower(argN)

                    case "detrend"

                        op.signal = cellfun(@(x) gen.make_row(squeeze(ep.detrend_(x, ndims(op.signal)))), op.signal, 'UniformOutput', false);

                    case "baseline"



                    case "subset"
                    
                    case "average"

                        subargN = varargin{iOper};
                        iOper = iOper + 1;
                        switch lower(subargN)

                            case 'trials'

                                avg_dim = 1;

                            case 'channels'

                                avg_dim = 2;

                            case 'all'

                                avg_dim = [1,2];

                        end

                        op.signal = cellfun(@(x) gen.make_row(squeeze(mean(x, avg_dim))), op.signal, 'UniformOutput', false);

                    case "pad"

                        subargN = varargin{iOper};
                        iOper = iOper + 1;
                        ep.isPadded_ = true;

                        switch lower(subargN)

                            case {'nan', 'nans'}

                                ep.pad_filler_ = NaN;

                            case {'zero', 'zeros'}

                                ep.pad_filler_ = 0;

                            otherwise

                                error('Epoch signal can only be padded with "zero" or "nan".');

                        end

                        ep.signal_length_after_padding_ = varargin{iOper};
                        iOper = iOper + 1;
                        assert(isnumeric(end_length), 'Pad length must be numeric and given as the 3rd input following the input "pad".');
                        
                    case "window"
                    case "fft"

                        sfreq = varargin{iOper};
                        iOper = iOper + 1;
                        if ~isnumeric(sfreq) && strcmp(sfreq, "sfreq")
                            sfreq = varargin{iOper};
                            iOper = iOper + 1;
                        end

                        if ep.isPadded_

                            tmp_signal = cellfun(@(x) paddata(x, ep.signal_length_after_padding_, Side='Both', FillValue=ep.pad_filler_, Dimension=ndims(x)), ...
                                op.signal, 'UniformOutput', false);

                        else

                            tmp_signal = op.signal;

                        end

                        for i = 1:height(op)

                            [op.amplitude{i}, op.phase{i}, fqN] = do_fft(op.signal{i}, ndims(tmp_signal{i}), sfreq);

                        end

                        ep.frequencies = fqN;

                    case "psd"

                        sfreq = varargin{iOper};
                        iOper = iOper + 1;
                        if ~isnumeric(sfreq) && strcmp(sfreq, "sfreq")
                            sfreq = varargin{iOper};
                            iOper = iOper + 1;
                        end

                        if nargin > iOper
                            opts = varargin{iOper};
                            if strcmp(opts, "options")

                                opts = varargin{iOper+1};
                                iOper = iOper + 2;
                            end

                        else
                            opts = {};
                        end

                        if ep.isPadded_

                            tmp_signal = cellfun(@(x) paddata(x, ep.signal_length_after_padding_, Side='Both', FillValue=ep.pad_filler_, Dimension=ndims(x)), ...
                                op.signal, 'UniformOutput', false);

                        else

                            tmp_signal = op.signal;

                        end

                        for i = 1:height(op)

                           [op.power{i}, fqN] = ep.do_psd(tmp_signal{i}, sfreq, ndims(tmp_signal{i}), opts{:});

                        end

                        ep.frequencies = fqN;

                    case "snr"

                        by_var = varargin{iOper};                        
                        if ~ismember(by_var, ["amplitude", "power"])
                            by_var = "amplitude";
                            if ~ismember(by_var, op.Properties.VariableNames)
                                by_var = "power";
                            end
                        else
                            iOper = iOper + 1;
                        end

                        subargN = varargin{iOper};
                        iOper = iOper + 1;
                        switch lower(subargN)
                            case {"neighbors"}

                                n_bins = varargin{iOper};
                                iOper = iOper + 1;

                                op.snr = cellfun(@(x) ...
                                    gen.apply_func_along_dimension(x, ndims(x), @ep.divide_by_n_neighbor_, n_bins), ...
                                    op.(by_var), 'UniformOutput', false);

                            case {"medfilt", "median"}

                                n_bins = varargin{iOper};
                                iOper = iOper + 1;

                                op.snr = cellfun(@(x) ...
                                    gen.apply_func_along_dimension(x, ndims(x), @ep.divide_by_medfilt_, n_bins), ...
                                    op.(by_var), 'UniformOutput', false);

                            case {"average"}
                            otherwise
                        end
                end

            end

            ep.data = op;
            ep.transform_steps_ = [ep.transform_steps_{:}, varargin{:}];
        end

        function win = get.epoch_win(ep)

            win = ep.EpochParameter.epoch_win;

        end

        function dt = get.dt(ep)

            dt = ep.Epoch.dt;

        end

        function t = get.timepoints(ep)

            win = ep.epoch_win;

            t = win(1):ep.dt:win(2);

            t = t(gen.ifwithin(t, ep.time_win));

        end

        function plot(ep, x, y, varargin)

            datN = ep.data;

            x = ep.(x);

            separate_figs_by = intersect(lower(datN.Properties.VariableNames), ep.separate_figure_by_levels_);
            n_figs_by_dim = varfun(@(x) length(unique(x)), datN(:,separate_figs_by), 'OutputFormat','uniform');
            n_figs = prod(n_figs_by_dim);

            n_subplot = length(unique(datN.name));

            for iFig = 1:n_figs

                figure;
                tileN = tiledlayout(n_subplot, 1, 'TileSpacing', 'none', 'Padding', 'none');

                for iSp = 1:n_subplot

                    axN = nexttile(tileN);
                    if ~isempty(datN.(y){iSp})
                        n_rows = size(datN.(y){iSp}, 1);
                        for j = 1:n_rows
                            plot(axN, x, squeeze(datN.(y){iSp}(j, :)), 'LineWidth', 1.5);
                            hold(axN, 'on');
                        end
                        hold(axN, 'off');
                    else
                        plot(axN, NaN, NaN);
                    end


                end

            end
        end
    end

    methods (Access = protected)

        function ep = retrieve(ep)

            eTbl = ep.Epoch;
            primary_keys = fetchtable(eTbl);
            uniq_keys = varfun(@(x) unique(x), primary_keys, 'OutputFormat', 'cell');
            n_uniq_keys = cellfun(@(x) length(unique(x)), uniq_keys);
            isDimExpand = n_uniq_keys ~= 1;

            op = table('Size', [height(primary_keys), sum(isDimExpand) + 1], ...
                'VariableTypes', [primary_keys.Properties.VariableTypes(isDimExpand), "cell"], ...
                'VariableNames', [primary_keys.Properties.VariableNames(isDimExpand), "signal"]);

            n_rows = height(primary_keys);

            if sum(isDimExpand)

                op(:, 1:sum(isDimExpand)) = primary_keys(:, isDimExpand);

            end

            for iKey = 1:n_rows

                qryN = table2struct(primary_keys(iKey,:));

                % Select trials to plot
                isValidN = ~isinf(getfield(fetch(ns.Epoch & qryN,'event_onset'),'event_onset')); % inf indicates that the trial was not shown
                trialsN = getfield(fetch(ns.DimensionCondition & ep.Epoch & qryN, 'trials'), 'trials'); % the actual trial numbers

                if ~isempty(ep.channels)
                    channel_qryN = sprintf("channel in (%s)", join(string(ep.channels),","));
                else 
                    channel_qryN = '';
                end

                % fetch data

                t_fetch = gen.Timer().start("Fetching epochs: %d of % d\n", iKey, n_rows);
                if isempty(channel_qryN)
                    datN = fetch(ns.EpochChannel & ep & qryN, 'signal');
                else
                    datN = fetch(ns.EpochChannel & ep & qryN & channel_qryN, 'signal');
                end

                t_fetch.stop("\t Fetching is complete. ");
                t_fetch.report();

                if ~isempty(ep.trials)

                    trials2fetch = intersect(ep.trials{iKey}, trialsN(isValidN));
                    isFetchTrials = ismember(trialsN(isValidN),trials2fetch);

                else

                    trials2fetch = trialsN(isValidN);
                    isFetchTrials = ones(1,sum(isValidN))==1;

                end


                datN = cat(3, datN(:).signal); % trial by timepoint by channel
                datN = permute(datN,[1, 3, 2]); % trial by channel by timepoint
                datN = datN(isFetchTrials,:,:);

                isT = ep.timepoints >= ep.time_win(1) & ep.timepoints <= ep.time_win(2);
                
                op.trials{iKey} = trials2fetch;
                op.signal{iKey} = datN(:,:,isT);

            end

            ep.data = op;

        end     

    end

    methods (Access = protected, Static)

        function [pow, frequencies] = do_psd(dat, fs, iDim, varargin)

                dim_size = size(dat);
                dim_size(iDim) = 1;
                slice_idx = arrayfun(@(x) 1:x, dim_size, 'UniformOutput',false);
                slice_idx{iDim} = ':';
                slice_idx = table2cell(combinations(slice_idx{:}));

                [pow, frequencies] = gen.apply_func_along_dimension(dat, iDim, @pspectrum, fs, varargin{:});

                select_idx = repmat({1},1, ndims(pow));
                select_idx{iDim} = ':';
                frequencies = squeeze(frequencies(select_idx{:}))';


            end

            function [amplitude, phase, frequencies] = do_fft(data, n, fs)
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

                amplitude = abs(fftResult(idx{:}))/N;
                phase = angle(fftResult(idx{:}));


        end

        function A = detrend_(A, dim)

            A = gen.apply_func_along_dimension(A,ndims(A),@detrend);
            return

            % A: Input array

            % Get the number of dimensions
            n_dims = ndims(A);
            if n_dims == 2, A = detrend(A); return; end

            % Create the permutation order
            perm_order = 1:n_dims;
            perm_order([dim, n_dims]) = perm_order([n_dims, dim]);

            % Permute the array
            A = permute(A, perm_order);

            % Detrend along the last dimension
            A = detrend(A);

            % Permute back to the original order
            A = permute(A, perm_order);

        end

        function denoised_vec = divide_by_n_neighbor_(vec, n)

            len = length(vec);
            noise_vec = zeros(size(vec));

            for i = 1:len

                if i <= n || i > len - n

                    noise_vec(i) = NaN;

                else

                    noise_vec(i) = mean(vec([(i-n):(i-1), (i+1):(i+n)]));

                end

            end

            denoised_vec = vec./noise_vec;

        end

        function denoised_vec = divide_by_medfilt_(vec, n)

            denoised_vec = vec ./ medfilt1(vec, n, 'truncate');

        end

    end
    
end