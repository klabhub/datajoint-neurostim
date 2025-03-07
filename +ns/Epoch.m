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

            signal = fetch(ns.CChannel & exp_tbl, 'signal');
            channels = fetch(ns.CChannel & exp_tbl, 'channel');

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

        function varargout = plot(eTbl, pv)

            % In the future add another method called get_epochs with
            % averaging options, basically it should include the code around after
            % '%fetch data' 

            arguments

                eTbl (1,1) ns.Epoch

                pv.channel = [] % if empty all channels
                pv.trials = [] % if empty, all trials
                pv.average = false  % if true,
                pv.collapse = 'trials' % can be 'channels', 'trials', 'all'

                pv.across = 'time' % can be 'frequency'
                pv.type = 'line' % can be 'raster'

                pv.time_window = [];
                pv.frequency_window = [];

                pv.frequency_output = 'amplitude' % can also be 'phase' 'power', 'decibel', 'snr', if the final, need to select snr method

                pv.snr_method = 'extract_neighbor' % can also be 'medfilt', 'average', or a custom function handle

            end

            % Decide on figures and subplotting
            primary_keys = fetchtable(eTbl);

            ctags = unique(primary_keys.ctag);
            n_ctag = length(ctags);
            
            etags = unique(primary_keys.etag);
            n_etag = length(etags);

            subs = unique(primary_keys.subject);
            n_sub = length(subs);  % each subject will have a different figure

            dims = unique(primary_keys.dimension);
            n_dim = length(dims); % each dimension will have a different figure

            conds = unique(primary_keys.name);
            n_cond = length(conds); % each level will be ploted at a different row in the same figure

            [~, dt] = sampleTime(ns.C & eTbl);
            epoch_win = getfield(getfield(fetch(ns.EpochParameter & eTbl,'pv'),'pv'),'epoch_win');
            t = epoch_win(1):dt:epoch_win(2);

            for iC = 1:n_ctag

                cN = ctags{iC};

                for iE = 1:n_etag

                    eN = etags{iE};
                
                    for iSub = 1:n_sub

                        subN = subs{iSub};

                        for iDim = 1:n_dim

                            dimN = dims{iDim};
                            iFig = sub2ind([n_ctag, n_etag, n_sub,n_dim], iC, iE, iSub, iDim); % figure number

                            % Conditions are subplotted
                            for iCond = 1:n_cond

                                condN = conds{iCond};
                                % fetch data

                                dat_qry = ns.Epoch.query_restricter( ...
                                    'ctag', cN, 'etag', eN, 'subject', subN, 'dimension', dimN, 'name', condN, ...
                                    'channel', pv.channel);
                                dat = fetch(ns.EpochChannel & eTbl & dat_qry, 'signal');
                                
                                % Select trials to plot                                
                                isValid = ~isinf(getfield(fetch(ns.Epoch & dat_qry,'event_onset'),'event_onset')); % inf indicates that the trial was not shown
                                trials = getfield(fetch(ns.DimensionCondition & eTbl & dat_qry, 'trials'), 'trials'); % the actual trial numbers

                                if ~isempty(pv.trials)

                                    trials2plot = intersect(pv.trials, trials(isValid));
                                    isPlotTrials = ismember(trials(isValid),trials2plot);

                                else

                                    isPlotTrials = ones(1,sum(isValid));

                                end

                                dat = cat(3, dat(:).signal); % trial by timepoint by channel
                                dat = permute(dat,[1, 3, 2]); % trial by channel by timepoint
                                dat = dat(isPlotTrials,:,:);

                                % Collapse
                                switch pv.collapse

                                    case 'channels'

                                        dat = squeeze(mean(dat,2));

                                    case 'trials'

                                        dat = squeeze(mean(dat,1));

                                    case 'all'

                                        dat = squeeze(mean(dat,[1,2]))';

                                    otherwise

                                        error("Only collapsing options are 'channels', 'trials', 'all'");

                                end
                                
                                % plots across time or frequency, the latter req
                                % uires more options

                                figure(iFig);
                                subplot(n_cond, 1, iCond);
                                switch pv.across

                                    case 'time'

                                        if isempty(pv.time_window)

                                            x = t;
                                            isInEpoch = ones(size(t));
                                        else

                                            isInEpoch = t>= pv.time_window(1) & t<=pv.time_window(2);
                                            x = t(isInEpoch);
                                        end
                                            x_unit = "ms";
                                            y_unit = "\muV";
                                            
                                    case 'frequency'

                                    otherwise

                                        error('Only can plot across time or frequency');

                                end

                                switch pv.type

                                    case 'line'

                                    case 'raster'

                                end

                                

                            end
    
                        end
                    end
                    
                end

            end



        %%

        end

    end


    methods (Static, Access = public)

        function s = query_restricter(varargin)

            if nargin==0 || rem(nargin,2)
                error("restricter takes in an even number of inputs such that \n \t restricter (field1, value1, ...).")
            end

            n_fields = nargin/2;

            for iField = 1:2:nargin

                if isempty(varargin{iField+1})

                    varargin(iField:(iField+1)) = [];

                end

            end

            s = struct(varargin{:});

        end

        function line_plot_(dat)
        end

        function raster_plot()
        end
       

    end

end