classdef RetrievedEpochs < matlab.mixin.Copyable
%% DESCRIPTION:
%
% Represents a collection of epochs (time segments of neural data) retrieved
% from an ns.Epoch data source. This class facilitates common operations
% like data transformation (detrending, baselining, FFT, PSD, SNR),
% subsetting (by time, channels, trials, frequency), and plotting.
%
% The core idea is to fetch raw epoch data based on specific criteria
% (channels, trials, time window) and store it along with relevant metadata
% in a table (`data`). Subsequent transformations modify this table,
% potentially adding new data columns (e.g., 'amplitude', 'power').
%
%% CONSTRAINTS:
% - The initial retrieval is based on a single ns.Epoch entry, implying
%   epochs share common parameters like sampling rate and original time
%   window, as asserted in the constructor.
% - Data within each cell of the `data` table (e.g., `data.signal{i}`) is
%   expected to have consistent dimensions initially, often
%   (trials x channels x time). Averaging operations can collapse
%   dimensions.
%
%% EXAMPLE
%
% ep = RetrievedEpochs(eTbl, 'channels', [*list of channels*], 'trials',
% {[*list of trial numbers per row in eTbl*],[...],...}, 'time_window', [*
% (1 x 2) range of timepoints to be extracted.*])
%
%% Inputs (Constructor)
%   eTbl (ns.Epoch)     : A scalar ns.Epoch object specifying the data source.
%                         Must correspond to a single unique key combination
%                         (e.g., one etag).
%   pv.channels (double): Optional. Vector of channel indices/IDs to retrieve.
%                         If empty or omitted, all channels are retrieved.
%   pv.trials (cell)    : Optional. Cell array specifying trials to retrieve
%                         for each potential row defined by `eTbl`'s primary keys.
%                         If empty or omitted, all valid trials are retrieved.
%                         Since eTbl has 1 row, this should be a 1x1 cell array.
%   pv.time_window (double): Optional. 1x2 vector [start_time, end_time]
%                         specifying the time window within the epoch to retrieve.
%                         If empty or omitted, uses the full epoch_win from
%                         ns.EpochParameter.
%
%% Properties (Public, Protected, Dependent)
% Public:
%   Epoch (ns.Epoch)    : Handle to the source ns.Epoch object.
%   EpochParameter      : Handle to the associated ns.EpochParameter object.
%
% Protected:
%   data (table)        : Core data container. Rows often represent different
%                         experimental conditions or groupings. Columns include
%                         metadata (from ns.Epoch primary keys) and cell arrays
%                         holding the numerical data ('signal', 'amplitude', etc.).
%   channels (double)   : Vector of channel IDs included in the data.
%   time_window (double): [start, end] time relevant to the current data (seconds).
%   frequency_window (double): [start, end] frequency window (Hz) after subsetting.
%   frequencies (double): Vector of frequency values (Hz) after FFT/PSD.
%   trials (cell)       : {1 x n_rows} cell array. Trial numbers for each row.
% Dependent:
%   timepoints (double) : Vector of time points for the data's time axis (seconds).
%   epoch_window (double): Original [start, end] epoch window from EpochParameter (seconds).
%   dt (double)         : Time step / sampling interval (seconds).
%   n_rows (double)     : Number of rows in the `data` table.
%
%% Methods (Key Public with Examples)
%   RetrievedEpochs     : Constructor - Initializes object and retrieves data.
%       Example: 
%           myEpoch = ns.Epoch & 'subject="S1"' & 'etag="VisualStim"';
%           rp = ns.RetrievedEpochs(myEpoch, channels=[1, 5, 10], time_window=[-0.5, 1.5]);
%
%   transform           : Applies data processing steps (detrend, baseline, FFT, etc.).
%       Examples:
%           rp = rp.transform('detrend'); 
%           rp = rp.transform('baseline', [-0.2, 0]); % Baseline correct using -200ms to 0ms
%           rp = rp.transform('average', 'trials'); % Average across trials
%           rp = rp.transform('pad', 'zero', 2048, 'psd', 'sfreq', rp.Epoch.sfreq, 'options', {{'Method','welch'}}); % Pad, then Welch PSD
%           rp = rp.transform('fft', 'sfreq', 1000, 'snr', 'amplitude', 'neighbors', 5); % FFT, then SNR on amplitude using +/- 5 neighbors
%           rp = rp.transform('split_trials', @(t) mod(t,2), 'Parity', {'Even','Odd'}); % Split by trial parity
%   subset              : Selects a subset of data based on criteria.
%       Example: 
%           subset_criteria.time_window = [0, 1.0];
%           subset_criteria.channels = [1, 5];
%           subset_criteria.frequency_window = [10, 30]; % Only if FFT/PSD done
%           subset_criteria.trials = {[1:10]}; % Assuming 1 row, select trials 1-10
%           rp = rp.subset(subset_criteria);
%   insert              : Adds new columns to the `data` table.
%       Example:
%           labels = {'ConditionA'}; % Assuming rp has 1 row currently
%           rp = rp.insert('ConditionLabel', labels);
%   plot                : Visualizes the data (e.g., signal vs. time, power vs. frequency).
%       Examples:
%           rp.plot('signal'); % Plot signal vs time
%           rp.plot('power', 'xlim', [5, 50], 'ylim', 'scale_figure'); % Plot power (5-50Hz), scale figures identically
%           rp.plot('amplitude', 'collate', 'subject'); % Plot amplitude, arrange figures by subject into one window

    properties

        Epoch (1,1) ns.Epoch
        EpochParameter (1,1) ns.EpochParameter        

    end

    properties (SetAccess = protected)

        data table
        channels (1,:) double = []
        time_window double = []
        frequency_window = []
        frequencies
        trials = {}

    end

    properties (Dependent)

        timepoints
        epoch_window
        dt
        n_rows

    end

    properties (Hidden)

        isPadded_ = false
        pad_filler_ = 'zeros';
        signal_length_after_padding_
        transform_steps_ = {}
        dimensions_ = struct(trials=1, channels=2, time=3, frequency=3)

    end

    properties (Hidden, Constant)

        separate_figure_by_levels_ = ["ctag", "etag", "dimension", "subject"]
        data_columns_ = ["signal", "amplitude", "power", "snr"]

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
            ep.time_window = pv.time_window;
            ep.trials = pv.trials;
            ep.retrieve;


        end

        function ep = transform(ep, varargin)


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