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
            % abs_onsets = abs_onsets(isVld);

            % conds = get(exp_tbl, parms.plugin_name, 'prm',parms.dimension,'what','data','atTrialTime',0)';
            % if isempty(conds), return; end
            % conds = conds(isVld);

            [t, ~] = tbl.get_sampling_times(ns.C & key);

            t_export = Timer().start("Exporting data from ns.CChannel...\n");

            signal = fetch(ns.CChannel & exp_tbl, 'signal');
            channels = fetch(ns.CChannel & exp_tbl, 'channel');

            t_export.stop("\tExporting is complete.");
            t_export.duration;

            t_sgm = Timer().start("Now segmenting...\n");

            sts = ns.SegmentedTimeSeries(t, abs_onsets(trials), parms.pv.epoch_win, horzcat(signal(:).signal)');
            sts.make();
            if parms.pv.baseline
                sts.baseline(parms.pv.baseline_win);
            end            

            ep = sts.epochs;

            t_sgm.stop("\t\tSegmentation complete.");
            t_sgm.duration;

            t_sub = Timer().start("\tNow submitting\n");

            epoch_tpl = mergestruct(key, struct( ...
                    event_onset = rel_onsets...
                    ));
            
            n_channels = size(ep,2);

            cepoch_tpl = mergestruct(repmat(key,n_channels,1), ...
                arrayfun(@(x) struct(channel=x.channel),channels), ...
                arrayfun(@(iChannel) struct(signal = squeeze(ep(:,iChannel,:))), [channels.channel]));

            insert(ns.Epoch, epoch_tpl);
            chunkedInsert(ns.EpochChannel, cepoch_tpl)

            t_sgm.stop("\t\tSubmission is complete.\n");
           
        end

    end

    methods (Static)

        function [t,dt] = get_sampling_times(ctbl)
            % Determine time and time step of this ns.C entry.
            arguments
                ctbl (1,1) {mustHaveRows(ctbl,1)}
            end
            t= double(fetch1(ctbl ,'time'));
            if numel(t)==3
                t= linspace(t(1),t(2),t(3))';
            end
            dt = mode(diff(t));
        end

    end

end