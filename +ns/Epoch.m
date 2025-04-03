%{
# Preprocessed epoch data, the actual epochs are stored in ns.EpochChannel

-> ns.C
-> ns.DimensionCondition
-> ns.EpochParameter

---
event_onset : BLOB # clock time to which the trials are locked, t0
%}

classdef Epoch < dj.Computed

    properties

        retrieved table = table() % the most recent retrieved data (the output of 'retrieve')
        channels (1,:) double = []
        frequencies (1,:) double = []

    end

    properties (Access = protected, Constant)

        args_for_retrieve_ = ["trials", "channels", "time_window"]

    end

    properties (Dependent)
        
        dt

    end

    methods (Access = protected)

        function makeTuples(eTbl, key)

            exp_tbl = ns.Experiment & key;
            if count(exp_tbl) == 0
                return;
            end

            pTbl = ns.EpochParameter & key;
            dims = fetch(ns.DimensionCondition & key, '*');

            abs_onsets = get(exp_tbl, 'cic','prm','firstFrame','atTrialTime',inf,'what','clocktime');
            rel_onsets = get(exp_tbl, pTbl.plugin_name, 'prm', 'startTime', 'what', 'trialtime');

            if isempty(abs_onsets), return; end
            isVld = ~isinf(rel_onsets) & ~isnan(rel_onsets);

            trials = 1:length(rel_onsets);
            trials = intersect(dims.trials, trials(isVld));

            [t, ~] = sampleTime(ns.C & key); % ns.C method
            t_export = gen.Timer().start("Exporting data from ns.CChannel...\n");

            signal = fetch(ns.CChannel & (ns.C & key), 'signal');
            ch_list = fetch(ns.CChannel & (ns.C & key), 'channel');

            t_export.stop("\tExporting is complete.");
            t_export.report();

            t_sgm = gen.Timer().start("Now segmenting...\n");

            sts = ns.SegmentedTimeSeries(t, abs_onsets(trials) + rel_onsets(trials), pTbl.epoch_win, horzcat(signal(:).signal)');
            sts.make();
            if pTbl.baseline
                sts.baseline(pTbl.baseline_win);
            end

            if pTbl.detrend
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
                arrayfun(@(x) struct(channel=x.channel),ch_list), ...
                    arrayfun(@(iChannel) struct(signal = squeeze(ep(:,iChannel,:))), [ch_list.channel]') ...
                    );

            insert(eTbl, epoch_tpl);
            chunkedInsert(ns.EpochChannel, cepoch_tpl)

            t_sub.stop("\t\tSubmission is complete.\n");
            t_sub.report();

        end        

    end
    methods

        function dt = get.dt(eTbl)

            c = fetch(ns.C & eTbl, 'time');
            dt = c(1).time;
            if numel(dt)==3
                dt = (dt(2)-dt(1))/dt(3);
            else
                dt = mode(diff(dt));
            end

        end

        function ep = retrieve(eTbl, varargin)

            ep = ns.RetrievedEpochs(eTbl, varargin{:});

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

 

end