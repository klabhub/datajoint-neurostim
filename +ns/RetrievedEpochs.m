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
        EpochChannel (1,1) ns.EpochChannel
        EpochParameter (1,1) ns.EpochParameter 
        Artifact (1,1) ns.Artifact
        ArtifactChannel (1,1) ns.ArtifactChannel
        C (1,1) ns.C

    end

    properties (SetAccess = protected)

        data table
        channels
        trials
        time_window double = []
        frequency_window = []
        frequencies
        artifacts struct = struct(keep=true, discard=false, fillnan=false, interpolate=false)
        badEpochs struct
    end

    properties (Dependent)

        timepoints
        epoch_length
        epoch_window
        dt
        n_rows
        data_columns
        temporal_data_columns
        spectral_data_columns

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
        temporal_data_columns_ = ["signal"];
        spectral_data_columns_ = ["amplitude", "power", "snr"];
    
    end

    methods

        function ep = RetrievedEpochs(eTbl, pv)

            arguments

                eTbl (1,1) ns.Epoch

                pv.channels = [] % if empty all channels
                pv.trials = {} % if empty, all trials
                pv.time_window = [];
                pv.artifacts = []

            end

            ep.Epoch = eTbl;
            ep.EpochParameter = ns.EpochParameter & eTbl;
            ep.EpochChannel = ns.EpochChannel & eTbl;
            ep.C = ns.C & eTbl;

            ep.Artifact = ns.Artifact & eTbl;
            ep.ArtifactChannel = ns.ArtifactChannel & eTbl;

            if isempty(pv.time_window)

                pv.time_window = ep.EpochParameter.epoch_win;

            end

            assert(ep.EpochParameter.nrows == 1, ...
                "retrieve(ns.Epoch,...) is designed for compatible epoch structures corresponding to a single etag.")


            ep.channels = pv.channels;
            ep.time_window = pv.time_window;
            ep.trials = pv.trials;
            if ~isempty(pv.artifacts)
                ep.artifacts = pv.artifacts;

                if islogical(struct2array(ep.artifacts))

                    assert(sum(struct2array(ep.artifacts))==1, "Only one of the fields in the property 'artifacts' must be set to true.");
                end
            end
            ep.retrieve;


        end

        function ep = transform(ep, varargin)

            iOper = 1;
            while iOper <= nargin-1

                argN = varargin{iOper};
                iOper = iOper + 1;
                switch lower(argN)

                    case "rereference"

                        if strcmp(varargin{iOper}, 'options')
                            
                            ref_opts = varargin{iOper + 1};
                            iOper = iOper + 2;
                            
                        else
                            ref_opts = {};
                        end

                        c_info = getfield(fetch(ep.C, 'info'),'info');

                        if isfield(c_info.layout, 'ChannelLocations')
                            ref_opts = [ref_opts, 'channel_locations', c_info.layout.ChannelLocations];
                        end                                            

                        ep.data.signal = cellfun(@(x) ns.rereference(x, ref_opts{:}), ep.data.signal, 'UniformOutput',false);

                    case "detrend"

                        ep.data.signal = cellfun(@(x) gen.make_row(squeeze(ep.detrend_(x))), ep.data.signal, 'UniformOutput', false);

                    case "baseline"

                        restricter = ns.EpochRestricter(time_window = varargin{iOper});
                        iOper = iOper + 1;

                        baselineN = ep.copy().subset(restricter).apply('signal', @mean, ep.dimensions_.time, 'omitnan');

                        % ep.data.signal = 
                        ep.data.signal = gen.make_column(arrayfun(@(i) ep.data.signal{i} - baselineN.data.signal{i}, 1:ep.n_rows, 'UniformOutput', false));

                    case "subset"

                        restricter = varargin{iOper};
                        if ~isa(restricter, 'ns.EpochRestricter')

                            error("'subset' command needs to be followed by a (1x1) structure array containing some or all of the following fieldL 'trials', 'channels', 'time_window'")
                        else
                            iOper = iOper + 1;
                        end

                        ep = ep.subset(restricter);

                    case "average"

                        subargN = varargin{iOper};
                        iOper = iOper + 1;
                        switch lower(subargN)

                            case 'trials'

                                assert(~isnan(ep.dimensions_.trials), "Trials were already collapsed.")
                                avg_dim = ep.dimensions_.trials;
                                ep.dimensions_.trials = NaN;
                                if ~isnan(ep.dimensions_.channels)
                                    ep.dimensions_.channels = 1;
                                end
                                ep.dimensions_.time = 2;

                            case 'channels'

                                assert(~isnan(ep.dimensions_.channels), "Channels were already collapsed.")
                                avg_dim = ep.dimensions_.channels;
                                ep.dimensions_.channels = NaN;
                                ep.dimensions_.time = 2;

                            case 'all'

                                assert(~any(isnan([ep.dimensions_.channels, ep.dimensions_.trials])), "Channels and/or trials were already collapsed.")
                                avg_dim = [ep.dimensions_.trials, ep.dimensions_.channels];
                                ep.dimensions_.trials = NaN;
                                ep.dimensions_.channels = NaN;
                                ep.dimensions_.time = 2;

                        end

                        for colN = ep.data_columns

                            ep.data.(colN) = cellfun(@(x) gen.make_row(squeeze(mean(x, avg_dim, 'omitmissing'))), ep.data.(colN), 'UniformOutput', false);

                        end

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
                        assert(isnumeric(ep.signal_length_after_padding_), 'Pad length must be numeric and given as the 3rd input following the input "pad".');

                    case "window"
                    case "fft"
                        
                        sfreq = 1000/ep.dt;

                        if ep.isPadded_

                            tmp_signal = cellfun(@(x) paddata(x, ep.signal_length_after_padding_, Side='Both', FillValue=ep.pad_filler_, Dimension=ndims(x)), ...
                                ep.data.signal, 'UniformOutput', false);

                        else

                            tmp_signal = ep.data.signal;

                        end

                        for i = 1:ep.n_rows

                            [ep.data.amplitude{i}, ep.data.phase{i}, fqN] = ep.do_fft(tmp_signal{i}, sfreq, ndims(tmp_signal{i}));

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
                                ep.data.signal, 'UniformOutput', false);

                        else

                            tmp_signal = ep.data.signal;

                        end

                        for i = 1:ep.n_rows

                            [ep.data.power{i}, fqN] = ep.do_psd(tmp_signal{i}, sfreq, ndims(tmp_signal{i}), opts{:});

                        end

                        ep.frequencies = fqN;

                    case "snr"

                        by_var = varargin{iOper};
                        if ~ismember(by_var, ["amplitude", "power"])
                            by_var = "amplitude";
                            if ~ismember(by_var, ep.data.Properties.VariableNames)
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

                                ep.data.snr = cellfun(@(x) ...
                                    gen.apply_func_along_dimension(x, ndims(x), @ep.divide_by_n_neighbor_, n_bins), ...
                                    ep.data.(by_var), 'UniformOutput', false);

                            case {"medfilt", "median"}

                                n_bins = varargin{iOper};
                                iOper = iOper + 1;

                                ep.data.snr = cellfun(@(x) ...
                                    gen.apply_func_along_dimension(x, ndims(x), @ep.divide_by_medfilt_, n_bins), ...
                                    ep.data.(by_var), 'UniformOutput', false);

                            case {"average"}
                            otherwise
                        end

                    case "split_trials"


                        trial_splitter = varargin{iOper}; % either n_row x 1 cell of trials or a function that outputs indices or boolean indices to select trials per group
                        iOper = iOper + 1;
 
                        col_name = varargin{iOper};
                        iOper = iOper + 1;

                        level_names = varargin{iOper};
                        iOper = iOper + 1;

                        ep.data = ep.split_trials_(ep.data, trial_splitter, col_name, level_names);


                end

            end

            isFunc = cellfun(@(x) isa(x, 'function_handle'), varargin);
            for i = 1:length(isFunc)
                if ~any(isFunc)
                    break;
                elseif isFunc(i)

                    varargin{i} = func2str(varargin{i});
                    
                end

            end
            ep.transform_steps_ = cat(2, ep.transform_steps_, varargin);

        end

        function ep = subset(ep, restricter)

            % NULL, ep.channels to be fixed, use EpochRestricter instead

            arguments

                ep
                restricter (1,1) ns.EpochRestricter

            end


            ep = restricter.restrict(ep);            

        end

        function find_bad_epochs(ep, varargin)

            dat_mat = cat(1, ep.data.signal{:});
            n_trl_per_row = cellfun(@(x) length(x), ep.trials);

            ep.badEpochs = ns.detect_outlier_epochs(dat_mat, 1000/ep.dt, varargin{:});

        end

        function plot_bad_epochs(ep, varargin)

            dat_mat = cat(1, ep.data.signal{:});

            % Plot by category

            bad_cats = string(fieldnames(ep.badEpochs));

            bad_cats = bad_cats(startsWith(bad_cats, "badBy") & cellfun(@(x) ~isempty(x), struct2cell(ep.badEpochs)));

            n_cat = length(bad_cats);

            fig = figure();
            tiles = tiledlayout(n_cat, 1, 'TileSpacing', 'compact', Padding='compact');
            
            for iCat = 1:n_cat

                datN = squeeze(mean(dat_mat(ep.badEpochs.(bad_cats(iCat)), ep.badEpochs.parameters.criterion_channels,:),2, 'omitnan'));
                tileN = nexttile(tiles);
                plot(tileN, ep.timepoints, datN);

                title(tileN, sprintf("%s (%d Epochs)", bad_cats(iCat), length(ep.badEpochs.(bad_cats(iCat)))));

            end
            title(tiles, "Bad Epochs By Category")

            % All epochs
            fig = figure();
            tiles = tiledlayout(2, 1, 'TileSpacing', 'compact', Padding='compact');
            datN = squeeze(mean(dat_mat(:, ep.badEpochs.parameters.criterion_channels,:),2, 'omitnan'));
            tileN = nexttile(tiles);
            plot(tileN, ep.timepoints, datN);
            title(tileN, sprintf("All Epochs (%d Epochs)", size(dat_mat,1)));

            % Good Epochs
            datN = squeeze(mean(dat_mat(setdiff(1:size(dat_mat,1), ep.badEpochs.all), ep.badEpochs.parameters.criterion_channels,:),2, 'omitnan'));
            tileN = nexttile(tiles);
            plot(tileN, ep.timepoints, datN);
            title(tileN, sprintf("Good Epochs (%d Epochs)", size(dat_mat,1)- numel(ep.badEpochs.all)));


        end
        
        function ep = select_rows(ep, idx)

            ep.data = ep.data(idx, :);

        end

        function ep = insert(ep, varargin)

            n_arg = nargin-1;

            for whArg = 1:2:n_arg

                valN = varargin{whArg+1};
                if isscalar(valN)
                    valN = repmat(valN,[ep.n_rows,1]);
                end
                ep.data.(varargin{whArg}) = valN;

            end

        end

        function ep = drop(ep, varargin)

            % drop data column

            n_arg = nargin-1;            

            for whArg = 1:n_arg
                
                ep.data.(varargin{whArg}) = [];

            end

        end

        function ep = rename(ep, old_names, new_names)

            ep.data = renamevars(ep.data, old_names, new_names);

        end

        function ep = replace(ep, varargin)

            vars = string(varargin{1:2:end});
            vals = varargin{2:2:end};
            isVarExist = ismember(vars, ep.data_columns);

            assert(all(isVarExist), "Data variable(s) %s do(es) not exist in the data table.", join(vars(isVarExist), ", "));

            for iVar = 1:length(vars)

                ep.data.(vars(iVar)) = vals(:,iVar);

            end

        end

        function n = get.n_rows(ep)

            n = height(ep.data);

        end

        function win = get.epoch_window(ep)

            win = ep.EpochParameter.epoch_win;

        end

        function dt = get.dt(ep)

            dt = ep.Epoch.dt;

        end

        function t = get.timepoints(ep)

            win = ep.epoch_window;

            t = win(1):ep.dt:win(2);

            t = t(gen.ifwithin(t, ep.time_window));

        end

        function l = get.epoch_length(ep)

            l = length(ep.timepoints);

        end

        function names = get.data_columns(ep)

            names = intersect(ep.data.Properties.VariableNames, ep.data_columns_);

        end

        function names = get.temporal_data_columns(ep)

            names = intersect(ep.data.Properties.VariableNames, ep.temporal_data_columns_);

        end

        function names = get.spectral_data_columns(ep)

            names = intersect(ep.data.Properties.VariableNames, ep.spectral_data_columns_);

        end

        function plot(ep, y, varargin)


            x_lim = false;
            y_lim = "scale_plot";
            c_lim = "scale_plot";
            plot_type = "line";
            isCollate = false;
            cmap = "copper";

            n_arg = nargin - 3;

            whArg = 1;
            while whArg <= n_arg

                argN = varargin{whArg};
                whArg = whArg + 1;

                switch argN

                    case "xlim"

                        subargN = varargin{whArg};
                        whArg = whArg + 1;

                        if isnumeric(subargN)

                            x_lim = subargN;

                        end

                    case "ylim"

                        subargN = varargin{whArg};
                        whArg = whArg + 1;

                        y_lim = subargN;

                    case "clim"

                        subargN = varargin{whArg};
                        whArg = whArg + 1;

                        c_lim = subargN;

                    case "collate"

                        collate_by = string(varargin{whArg}); %by which variable
                        n_coll_var = length(collate_by);
                        whArg = whArg + 1;
                        assert(n_coll_var<=2, "Only 1 or 2 variables/columns are collateable.");
                        isCollate = true;

                    case "type"

                        plot_type = varargin{whArg};
                        whArg = whArg + 1;

                    case "colormap"

                        cmap = varargin{whArg};
                        whArg = whArg + 1;                  

                end

            end

            if strcmp(y, "signal")

                x_label = "Time (milliseconds)";
                x = "timepoints";
                y_label = "$\muV$";

            else

                x = "frequencies";
                x_label = "Frequency (Hz)";

            end

            if strcmp(y, "amplitude")
                y_label = "Amplitude (\muV)";

            elseif strcmp(y, "phase")
                y_label = "Phase angle (radians)";
            elseif strcmp(y, "power")
                y_label = "Power (\muV^2)";
            elseif strcmp(y, "snr")
                y_label = "SNR";
            end

            datN = ep.data;

            x = ep.(x);

            separate_figs_by = intersect(lower(datN.Properties.VariableNames), ep.separate_figure_by_levels_);
            n_figs_by_dim = varfun(@(x) length(unique(x)), datN(:,separate_figs_by), 'OutputFormat','uniform');
            if isempty(n_figs_by_dim); n_figs_by_dim = 1; end
            n_figs = prod(n_figs_by_dim);

            n_subplot = length(unique(datN.name));

            ax = squeeze(cell([n_figs_by_dim, n_subplot, 1]));
            figs = squeeze(cell([n_figs_by_dim, 1]));
            fig_subs = arrayfun(@(x) 1:x, n_figs_by_dim, 'UniformOutput', false);
            fig_subs = table2cell(combinations(fig_subs{:}));
            % fig_subsN  = squeeze(cell(size(n_figs_by_dim)));
            if isscalar(n_figs_by_dim), n_figs_by_dim = [n_figs_by_dim, 1]; end

            for iFig = 1:n_figs

                % [fig_subsN{:}] = ind2sub(n_figs_by_dim, iFig);
                fig_subsN = fig_subs(iFig,:);
                figs{fig_subsN{:}} = figure('Visible', 'off');
                colormap(figs{fig_subsN{:}}, cmap);
                tileN = tiledlayout(figs{fig_subsN{:}},n_subplot, 1, 'TileSpacing', 'loose', 'Padding', 'compact');

                datN = ep.data((1:n_subplot) + (iFig-1)*n_subplot,:);
                for iSp = 1:n_subplot

                    axN = nexttile(tileN);


                    if ~isempty(datN.(y){iSp})

                        if strcmp(plot_type, "line")
                            n_row = size(datN.(y){iSp}, 1);
                            for j = 1:n_row
                                plot(axN, x, squeeze(datN.(y){iSp}(j, :)), 'LineWidth', 1.5);
                                hold(axN, 'on');
                            end
                        elseif strcmp(plot_type, "raster")

                            if isnan(ep.dimensions_.channels)
                                
                                y_data = datN.trials{iSp};

                            elseif isnan(ep.dimensions_.trials)

                                y_data = datN.channels;

                            end

                            imagesc(axN, 'XData', x, 'YData', y_data, 'CData', squeeze(datN.(y){iSp}));

                            ylim(gen.range(y_data));
                            
                            if strcmp(c_lim, "scale_plot")
                                maxima = max(squeeze(datN.(y){iSp}),[],2);
                                maxima = maxima(~isnan(maxima));
                                maxima_scaled = gen.robust_z(maxima, [], "std");

                                clim(axN, [0, max(maxima(maxima_scaled < 1.96))])

                            elseif isnumeric(c_lim)

                                clim(axN,c_lim);

                            end

                        end

                        if isnumeric(x_lim)

                            xlim(axN, x_lim);

                        end


                        hold(axN, 'off');
                    else
                        plot(axN, NaN, NaN);
                    end
                    if iSp < n_subplot, set(axN, 'XTick',[], 'XTickLabel',[]); end
                    title(axN, datN.name{iSp}, 'Interpreter', 'latex');

                    ax{fig_subsN{:}, iSp} = axN;
                end

                if strcmp(y_lim, "scale_figure")

                    y_scale = vertcat(datN.(y){:});
                    if isnumeric(x_lim)
                        y_scale = y_scale(:, gen.ifwithin(x,x_lim));
                    end
                    y_limN = gen.range(y_scale,'all');
                    cellfun(@(axN) ylim(axN, y_limN), ax(fig_subsN{:},:));

                end

                if isnumeric(y_lim)

                    cellfun(@(axN) ylim(axN, y_lim), ax(fig_subsN{:},:));

                end


                ylabel(tileN, y_label);
                xlabel(tileN, x_label, 'Interpreter', 'latex');
                if ~isempty(separate_figs_by)
                    fig_titleN = join(arrayfun(@(col) sprintf("%s: %s", col, datN.(col){1}), separate_figs_by),",   ");
                    title(tileN, fig_titleN);
                end

            end

            if isCollate

                isCollateDim = ismember(separate_figs_by, collate_by);
                n_figs_to_collate = prod(n_figs_by_dim(isCollateDim));

                n_fig_dims = length(n_figs_by_dim);
                collate_indices = repmat({':'},1,n_fig_dims);
                for iDim = 1:n_fig_dims

                    if ~isCollateDim(iDim)

                        collate_indices{iDim} = 1:n_figs_by_dim(iDim);

                    end

                end
                collate_indices = table2cell(combinations(collate_indices{:}));

                for iFig = 1:size(collate_indices,1)

                    figsN = squeeze(figs(collate_indices{iFig,:}));
                    axN = squeeze(ax(collate_indices{iFig,:},:));

                    if n_coll_var == 1

                        sp_size = [round(sqrt(n_figs_to_collate)), ceil(sqrt(n_figs_to_collate))];
                        new_figs = cell(sp_size);

                        new_figs(1:numel(figsN)) = figsN(:);

                        axN = squeeze(permute(reshape(permute(cat(1,axN, cell(prod(sp_size)-size(axN,1),n_subplot)), [2:ndims(axN), 1]), [n_subplot, sp_size]), [ndims(axN), ndims(axN)+1, 1:(ndims(axN)-1)]));
                        figsN = squeeze(new_figs);



                    elseif n_coll_var == 2

                        sp_dim_order = arrayfun(@(x) find(ismember(separate_figs_by(isCollateDim), x)), collate_by);
                        % figsN = permute(figsN, sp_dim_order);
                        % axN = permute(axN,[sp_dim_order, 3]);


                    end

                    newFig = gen.combine_figures_into_subplots(figsN, axN);
                end
            else

                cellfun(@(figN) set(figN, 'Visible','on'), figs);


            end
            close
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

            n_row = height(primary_keys);

            if sum(isDimExpand)

                op(:, 1:sum(isDimExpand)) = primary_keys(:, isDimExpand);

            end

            for iKey = 1:n_row

                qryN = table2struct(primary_keys(iKey,:));

                isValidN = ~isinf(getfield(fetch(ns.Epoch & qryN,'event_onset'),'event_onset')); % inf indicates that the trial was not shown
                trialsN = getfield(fetch(ns.DimensionCondition & ep.Epoch & qryN, 'trials'), 'trials'); % the actual trial numbers

                if isempty(ep.channels) || (iscell(ep.channels) && isempty(ep.channels{iKey}))
                    
                    channel_qryN = '';
                
                else

                    if iscell(ep.channels), ch_list = ep.channels{iKey};
                    else, ch_list = ep.channels;
                    end
                    channel_qryN = sprintf("channel in (%s)", join(string(ch_list),","));
                    
                end

                % fetch data

                t_fetch = gen.Timer().start("Fetching epochs: %d of % d\n", iKey, n_row);
                artN = fetch(ep.Artifact & qryN);
                if isempty(channel_qryN)

                    ep_ch_qryN = ep.EpochChannel & qryN;
                    art_ch_qryN = ep.ArtifactChannel & qryN;
                    
                else

                    ep_ch_qryN = ep.EpochChannel & qryN & channel_qryN;
                    art_ch_qryN = ep.ArtifactChannel & qryN & channel_qryN;                    
                
                end

                datN = fetch(ep_ch_qryN, 'signal');
                
                t_fetch.stop("\t Fetching is complete. ");
                t_fetch.report();

                if ~isempty(ep.trials)

                    trials_to_fetch = intersect(ep.trials{iKey}, trialsN(isValidN));
                    isFetchTrials = ismember(trialsN(isValidN),trials_to_fetch);

                else

                    trials_to_fetch = trialsN(isValidN);
                    isFetchTrials = ones(1,sum(isValidN))==1;

                end

                %% Include artifacts or not
                isFillNaN = ones([length(isFetchTrials), height(datN), ep.epoch_length])==0;

                for fldN = string(fieldnames(ep.artifacts))'


                    if isfield(ep.artifacts, fldN) && islogical(ep.artifacts.(fldN)) 
                        
                        if ep.artifacts.(fldN)

                            atagN = unique(string({artN.atag}));

                        else, continue;
                        end

                    elseif isstring(ep.artifacts.(fldN))

                        atagN = ep.artifacts.(fldN);                        

                    end

                    switch fldN

                        case 'keep'

                            continue;

                        case {'discard', 'fillnan'}


                            atag_qryN = sprintf('atag in ("%s")', join(atagN,'","'));
                            
                            %% Check trial specific artifacts across all channels
                            tbl_qryN = ns.Artifact & artN & atag_qryN;                            
                            art_tblN = fetchtable(tbl_qryN, 'trial');
                            isTrl2DiscardN = arrayfun(@(x) double(x), art_tblN.trial, 'UniformOutput', false);
                            isTrl2DiscardN = ismember(trials_to_fetch, [isTrl2DiscardN{:}]);
                            
                            if strcmp(fldN, 'discard')
                                isFetchTrials(isTrl2DiscardN) = false;
                            else
                                isFillNaN(isTrl2DiscardN,:,:) = true;
                            end

                            art_ch_qryNN = art_ch_qryN & atag_qryN;
                                                     
                            if count(art_ch_qryNN)

                                ch_tblN = fetchtable(art_ch_qryNN,'*');
                                for ii = 1:height(ch_tblN)
                                    
                                    isChN = ch_tblN.channel(ii) == [datN.channel];
                                    if ch_tblN.Properties.VariableTypes(strcmp(ch_tblN.Properties.VariableNames,'trial'))=="double"
                                        trials_w_artN = ch_tblN.trial(ii,:);
                                    else
                                        trials_w_artN = ch_tblN.trial{ii};
                                    end
                                    isTrl2DiscardN = ismember(trials_to_fetch, trials_w_artN);
                                    if isempty(trials_w_artN) % faulty channels
                                        isTrl2DiscardN = ones(size(trials_to_fetch))==1;
                                    end
                                    % add start stop
                                    isFillNaN(isTrl2DiscardN,isChN, :) = true;

                                end

                            end

                        case 'interpolate'
                        otherwise
                            error("err msg")
                    end

                end
                
                datN = cat(3, datN(:).signal); % trial by timepoint by channel
                datN = permute(datN,[1, 3, 2]); % trial by channel by timepoint
                datN(isFillNaN) = NaN;
                datN = datN(isFetchTrials,:,:);
                

                isT = ep.timepoints >= ep.time_window(1) & ep.timepoints <= ep.time_window(2);

                op.trials{iKey} = trials_to_fetch(isFetchTrials);
                op.channels{iKey} = arrayfun(@(s) s.channel, fetch(ep_ch_qryN,'channel'));
                op.signal{iKey} = datN(:,:,isT);

            end

            ep.data = op;
            ep.channels = op.channels;
            ep.trials = op.trials;

        end

        function ep = apply(ep, var, func, varargin)

            for i = 1:ep.n_rows

                ep.data.(var){i} = func(ep.data.(var){i}, varargin{:});

            end

        end

        function ep = apply_func2rows(ep, func, varargin)

            for i = 1:ep.n_rows

                ep(i,:) = func(ep(i,:), varargin{:});

            end

        end



    end
    
    methods (Static)

        function ep = merge(ep1, ep2, variable_name, values)

            ep = ep1.copy();
            ep.insert(variable_name, values(1));            

            ep.data = vertcat(ep.data, ep2.copy().insert(variable_name, values(2)).data);
            ep.separate_figure_by_levels_ = [ep.separate_figure_by_levels_, variable_name];
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

        function A = detrend_(A)

            % warning('off', 'MATLAB:detrend:PolyNotUnique');
            n_dim = ndims(A);
            A = permute(A, circshift(1:n_dim,1));
            A = detrend(A,1,'omitnan');
            A = permute(A, circshift(1:n_dim,-1));

            % A = gen.apply_func_along_dimension(A,ndims(A),@detrend, 'omitnan');
            % warning('on', 'MATLAB:detrend:PolyNotUnique');
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

        function new_data = split_trials_(data, trials, col_name, level_names)

            n_levels = numel(level_names);



            var_names = data.Properties.VariableNames;
            isDataColumn = ismember(var_names,ns.RetrievedEpochs.data_columns_);
            copy_vars = var_names(~isDataColumn);
            data_cols = var_names(isDataColumn);
            new_data = repelem(data(:,copy_vars),n_levels,1);
            new_data.(col_name) = repmat(gen.make_column(level_names),height(data),1);
            new_data(:,data_cols) = repmat(cell(1),height(new_data),1);

            idx_selector = repmat({':'}, 1, ndims(data.(data_cols{1}){1}));

            iRow = 1;
            for i = 1:height(data)

                for iLvl = 1:n_levels

                    trlN = new_data{iRow,"trials"}{1};
                    if isa(trials, "function_handle")

                        isTrl = trials(trlN);
                        trlN = isTrl == (iLvl-1);

                    elseif isa(trials, "cell")

                        trls = trials{iLvl};
                        if isa(trls, "function_handle")

                            trls = trls(trlN);

                        end

                        if islogical(trls)

                            trlN = trls;

                        else

                            trlN = ismember(trlN,trls);

                        end

                    end

                    idx_selector{1} = trlN;
                    new_data.trials{iRow} = new_data.trials{iRow}(trlN);

                    for iCol = 1:sum(isDataColumn)

                        new_data(iRow, data_cols{iCol}) = {data.(data_cols{iCol}){i}(idx_selector{:})};
                    end

                    iRow = iRow + 1;

                end

            end


        end
    end

end