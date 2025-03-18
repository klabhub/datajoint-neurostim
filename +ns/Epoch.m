%{
# Preprocessed epoch data, the actual epochs are stored in ns.EpochChannel

-> ns.C
-> ns.DimensionCondition
-> ns.EpochParameter

---
event_onset : BLOB # clock time to which the trials are locked, t0
%}

classdef Epoch < dj.Computed

    methods (Access = protected)

        function makeTuples(tbl, key)

            exp_tbl = ns.Experiment & key;
            if count(exp_tbl) == 0
                return;
            end
            parms = fetch(ns.EpochParameter & key, '*');
            dims = fetch(ns.DimensionCondition & key, '*');


            abs_onsets = get(exp_tbl, parms.plugin_name, 'prm', 'startTime', 'what', 'clocktime');
            rel_onsets = get(exp_tbl, parms.plugin_name, 'prm', 'startTime', 'what', 'trialtime');

            if isempty(abs_onsets), return; end
            isVld = ~isinf(abs_onsets) & ~isnan(abs_onsets);

            trials = 1:length(rel_onsets);
            trials = intersect(dims.trials, trials(isVld));

            [t, ~] = sampleTime(ns.C & key); % ns.C method

            t_export = gen.Timer().start("Exporting data from ns.CChannel...\n");

            signal = fetch(ns.CChannel & (ns.C & key), 'signal');
            channels = fetch(ns.CChannel & (ns.C & key), 'channel');

            t_export.stop("\tExporting is complete.");
            t_export.report();

            t_sgm = gen.Timer().start("Now segmenting...\n");

            sts = ns.SegmentedTimeSeries(t, abs_onsets(trials), parms.pv.epoch_win, horzcat(signal(:).signal)');
            sts.make();
            if parms.pv.baseline
                sts.baseline(parms.pv.baseline_win);
            end

            if parms.pv.detrend
                sts.detrend();
            end

            ep = sts.epochs;

            t_sgm.stop("\t\tSegmentation complete.");
            t_sgm.report();

            t_sub = gen.Timer().start("\tNow submitting to the server\n");

            epoch_tpl = mergestruct(key, struct( ...
                event_onset = rel_onsets(dims.trials)...
                ));

            n_channels = size(ep,2);

            cepoch_tpl = mergestruct(repmat(key,n_channels,1), ...
                arrayfun(@(x) struct(channel=x.channel),channels), ...
                arrayfun(@(iChannel) struct(signal = squeeze(ep(:,iChannel,:))), [channels.channel]));

            insert(ns.Epoch, epoch_tpl);
            chunkedInsert(ns.EpochChannel, cepoch_tpl)

            t_sub.stop("\t\tSubmission is complete.\n");
            t_sub.report();

        end

    end
    methods (Access = public)

        function [op, t] = retrieve(eTbl, pv)

            arguments

                eTbl (1,1) ns.Epoch

                pv.channel = [] % if empty all channels
                pv.trials = {} % if empty, all trials
                pv.time_window = [];


            end
           
            epoch_win = getfield(getfield(fetch(ns.EpochParameter & eTbl,'pv'),'pv'),'epoch_win');

            if isempty(pv.time_window)

                pv.time_window = epoch_win;

            end


            % Decide on figures and subplotting
            primary_keys = fetchtable(eTbl);
            uniq_keys = varfun(@(x) unique(x), primary_keys, 'OutputFormat', 'cell');
            n_uniq_keys = cellfun(@(x) length(unique(x)), uniq_keys);
            isDimExpand = n_uniq_keys ~= 1;
            uniq_keys = cell2struct(uniq_keys', primary_keys.Properties.VariableNames);

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
                isValid = ~isinf(getfield(fetch(ns.Epoch & qryN,'event_onset'),'event_onset')); % inf indicates that the trial was not shown
                trials = getfield(fetch(ns.DimensionCondition & eTbl & qryN, 'trials'), 'trials'); % the actual trial numbers

                if ~isempty(pv.channel)
                    channel_qryN = sprintf("channel in (%s)", join(string(pv.channel),","));
                else, channel_qryN = '';
                end

                % fetch data

                t_fetch = gen.Timer().start("Fetching epochs: %d of % d\n", iKey, n_rows);
                if isempty(channel_qryN)
                    datN = fetch(ns.EpochChannel & eTbl & qryN, 'signal');
                else
                    datN = fetch(ns.EpochChannel & eTbl & qryN & channel_qryN, 'signal');
                end

                t_fetch.stop("\t Fetching iscomplete. ");
                t_fetch.report();

                if ~isempty(pv.trials)

                    trials2fetch = intersect(pv.trials{iKey}, trials(isValid));
                    isFetchTrials = ismember(trials(isValid),trials2fetch);

                else

                    isFetchTrials = ones(1,sum(isValid))==1;

                end


                datN = cat(3, datN(:).signal); % trial by timepoint by channel
                datN = permute(datN,[1, 3, 2]); % trial by channel by timepoint
                datN = datN(isFetchTrials,:,:);

                t = linspace(epoch_win(1),epoch_win(2), size(datN,3));
                isT = t >= pv.time_window(1) & t <= pv.time_window(2);
                datN = datN(:,:,isT);
                op.signal{iKey} = datN;
                op.time{iKey} = t(isT);

            end          

        end

        function varargout = plot(eTbl, pv)

            % In the future add another method called get_epochs with
            % averaging options, basically it should include the code around after
            % '%fetch data'

            arguments

                eTbl (1,1) ns.Epoch

                pv.channel = [] % if empty all channels
                pv.trials = [] % if empty, all trials
                pv.average = 'trials' % can be 'channels', 'trials', 'all'

                pv.across = 'time' % can be 'frequency' or 'time-frequency'
                pv.type = 'line' % can be 'raster'

                pv.time_window = [];
                pv.frequency_window = [];

                pv.frequency_output = 'amplitude' % can also be 'phase' 'power', 'decibel', 'snr', if the final, need to select snr method

                pv.snr_method = 'extract_neighbor' % can also be 'medfilt', 'average', or a custom function handle

            end

            assert(ismember(pv.average, {'channels', 'trials', 'all'}), ...
                "ns.Epoch can only plot averages across 'channels', 'trials', or 'all'.");

            [datN, t] = retrieve(eTbl, "channel", pv.channel, "trials", pv.trials, "time_window", pv.time_window);
            datN = ns.Epoch.transform(datN, "average", pv.average);

            separate_figs_by = intersect(lower(datN.Properties.VariableNames), ["subject", "etag", "ctag", "dimension"]);
            n_figs_by_dim = varfun(@(x) length(unique(x)), datN(:,separate_figs_by), 'OutputFormat','uniform');
            n_figs = prod(n_figs_by_dim);

            n_subplot = length(unique(datN.name));

            for iFig = 1:n_figs

                figure;
                tileN = tiledlayout(n_subplot, 1, 'TileSpacing', 'none', 'Padding', 'none');

                for iSp = 1:n_subplot

                    axN = nexttile(tileN);
                    if ~isempty(datN.signal{iFig})
                        n_rows = size(datN.signal{iFig}, 1);
                        for j = 1:n_rows
                            plot(ax, t, datN.signal{iFig}(j, :), 'LineWidth', 1.5);
                            hold(ax, 'on');
                        end
                        hold(ax, 'off');
                    else
                        plot(ax, NaN, NaN);
                    end


                end


            end


            %%
            
        end

        function trials = get_trial_numbers(eTbl)

                isValid = fetch(eTbl,'event_onset'); % inf indicates that the trial was not shown
                n_key = length(isValid);
                trials = fetchtable(eTbl);
                for i = 1:n_key

                    trlN = fetch(ns.DimensionCondition & isValid(i), 'trials'); % the actual trial numbers
                    trials.trial{i} = [trlN.trials(~isinf(isValid(i).event_onset))];

                end

            end

    end


    methods (Static, Access = public)

        function op = transform(op, varargin)

            iOper = 1;

            while iOper <= nargin-1

                argN = varargin{iOper};
                iOper = iOper + 1;
                switch lower(argN)
                    
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

                        op.signal = cellfun(@(x) ns.Epoch.make_row_if_column_(squeeze(mean(x, avg_dim))), op.signal, 'UniformOutput', false);

                    case "pad"

                        subargN = varargin{iOper};
                        iOper = iOper + 1;

                        switch lower(subargN)

                            case {'nan', 'nans'}

                                filler = NaN;

                            case {'zero', 'zeros'}

                                filler = 0;

                            otherwise

                                error('Epoch signal can only be padded with "zero" or "nan".');

                        end

                        end_length = varargin{iOper};
                        iOper = iOper + 1;
                        assert(isnumeric(end_length), 'Pad length must be numeric and given as the 3rd input following the input "pad".');
                        op.signal = cellfun(@(x) paddata(x, end_length, Side='Both', FillValue=filler, Dimension=ndims(x)), ...
                            op.signal, 'UniformOutput', false);

                    case "window"
                    case "fft"

                        sfreq = varargin{iOper};
                        iOper = iOper + 1;
                        if ~isnumeric(sfreq) && strcmp(sfreq, "sfreq")
                            sfreq = varargin{iOper};
                            iOper = iOper + 1;
                        end

                        for i = 1:height(op)

                            [op.amplitude{i}, op.phase{i}, op.frequency{i}] = do_fft(op.signal{i}, ndims(op.signal{i}), sfreq);

                        end

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

                        for i = 1:height(op)

                            [op.power{i}, op.frequency{i}] = do_psd(op.signal{i}, sfreq, ndims(op.signal{i}), opts{:});

                        end

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
                                    gen.apply_func_along_dimension(x, ndims(x), @ns.Epoch.divide_by_n_neighbor_, n_bins), ...
                                    op.(by_var), 'UniformOutput', false);

                            case {"medfilt", "median"}

                                n_bins = varargin{iOper};
                                iOper = iOper + 1;

                                op.snr = cellfun(@(x) ...
                                    gen.apply_func_along_dimension(x, ndims(x), @ns.Epoch.divide_by_medfilt_, n_bins), ...
                                    op.(by_var), 'UniformOutput', false);

                            case {"average"}
                            otherwise
                        end
                end

            end

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

                amplitude = abs(fftResult(idx{:}));
                phase = angle(fftResult(idx{:}));


            end

        end

        function line_plot_(dat)
        end

        function raster_plot()
        end

    end

    methods (Static, Access = protected)

        function arr = make_row_if_column_(arr)

            if iscolumn(arr), arr = arr'; end
        end

        function arr = make_column_if_row_(arr)

            if isrow(arr), arr = arr'; end
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