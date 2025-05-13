%{
# Preprocessed epoch data, the actual epochs are stored in ns.EpochChannel

-> ns.C
-> ns.Dimension
-> ns.EpochParameter
---
event_onset : BLOB # clock time to which the trials are locked, t0
noisy_epochs : longblob # outlier epoch flags
%}

classdef Epoch < dj.Computed & dj.DJInstance

    properties (Dependent)

        C
        channels
        timepoints
        frequencies
        artifacts
        trials      
        conditions

        keySource

    end

    properties (Access = protected)

        trials_ = ns.Epoch_trials()
        conditions_

    end

    methods

        function v = get.keySource(~)

            v = ns.C * ns.Dimension * ( ns.Experiment *proj(ns.EpochParameter,'dimension'));
        
        end
    end

    methods (Access = protected)

        function makeTuples(eTbl, key)

            exp_tbl = ns.Experiment & key;           

            epTbl = ns.EpochParameter & key;
            cTbl = ns.C & key;
            dims = fetch(ns.DimensionCondition & key, '*');

            abs_onsets = get(exp_tbl, 'cic','prm','firstFrame','atTrialTime',inf,'what','clocktime');
            rel_onsets = get(exp_tbl, epTbl.plugin_name, 'prm', 'startTime', 'what', 'trialtime');
            
            if isempty(abs_onsets), return; end
            isVld = ~isinf(rel_onsets) & ~isnan(rel_onsets);

            trial_no = sort(cat(1,dims.trials));%1:length(rel_onsets);
            trial_no = trial_no(isVld);
            % trial_no = intersect(dims.trials, trial_no(isVld));
            onsets = abs_onsets(trial_no) + rel_onsets(trial_no);
            
            t = cTbl.timepoints;
            dt = cTbl.dt;
            % [t, dt] = sampleTime(ns.C & key); % ns.C method
            ep_timepoints = epTbl.epoch_win(1):dt:epTbl.epoch_win(2);

            t_export = gen.Timer().start("Exporting data from ns.CChannel...\n");

            signal = fetch(ns.CChannel & cTbl, 'signal');
            signal = horzcat(signal(:).signal)';
            ch_list = cTbl.channels;

            t_export.stop("\tExporting is complete.");
            t_export.report();

            t_sgm = gen.Timer().start("Now segmenting...\n");

            
            % sts = ns.SegmentedTimeSeries(t, abs_onsets(trials) + rel_onsets(trials), epTbl.epoch_win, horzcat(signal(:).signal)');
            % sts.make();
            % if epTbl.baseline
            %     sts.baseline(epTbl.baseline_win);
            % end
            % 
            % if epTbl.detrend
            %     sts.detrend();
            % end
            % 
            % ep = sts.epochs;

            ep = segment();
            pv = getfield(fetch(epTbl, 'pv'), 'pv');
            
            if isfield(pv, 'detrend')

                ep = eTbl.detrend(ep);

            end

            if isfield(pv, 'baseline') && pv.baseline

                ep = eTbl.baseline(ep, gen.ifwithin(ep_timepoints, pv.baseline_win));

            end

            if isfield(pv, 'rereference')

                ep(:, cTbl.artifacts.all,:) = NaN; % exclude noisy channels
                ep = eTbl.rereference(ep, cTbl.channel_coordinates, pv.ref_opts{:});

            end

            if isfield(pv, 'find_outliers') && pv.find_outliers

                if isfield(pv, 'outlier_opts')
                    outlier_opts = pv.outlier_opts;
                else
                    outlier_opts = {};
                end
                noisy_epochs = ns.detect_outlier_epochs(ep, 1000/cTbl.dt, outlier_opts{:});
            else
                noisy_epochs = [];
            end

            t_sgm.stop("\t\tSegmentation complete.");
            t_sgm.report();

            t_sub = gen.Timer().start("\tNow submitting to the server\n");

            epoch_tpl = mergestruct(key, struct( ...
                noisy_epochs = noisy_epochs,...
                event_onset = rel_onsets(trial_no)...
                ));

            n_channels = size(ep,2);

            % cepoch_tpl = mergestruct(repmat(key,n_channels,1), ...
            %     arrayfun(@(x) struct(channel=x.channel),ch_list), ...
            %         arrayfun(@(iChannel) struct(signal = squeeze(ep(:,iChannel,:))), [ch_list.channel]') ...
            %         );

            cepoch_tpl = mergestruct(repmat(key,n_channels,1), ...
                arrayfun(@(x) struct(channel=x),ch_list), ...
                    arrayfun(@(iChannel) struct(signal = squeeze(ep(:,iChannel,:))), ch_list) ...
                    );

            insert(eTbl, epoch_tpl);
            chunkedInsert(ns.EpochChannel, cepoch_tpl)

            t_sub.stop("\t\tSubmission is complete.\n");
            t_sub.report();

            function ep = segment()
                
                % Get indices for each epoch start, onset, and end
                
                nTrials = length(onsets);
                nChannels = size(signal,1);
                nEpochSamples = length(ep_timepoints);
                nSample = size(t, 2);
                iStart = arrayfun(@(i) gen.absargmin(t - (onsets(i) + epTbl.epoch_win(1))), 1:nTrials);
                % iOnset = iStart + find(ep_timepoints >= 0, 1, 'first') - 1;
                iEnd = iStart + nEpochSamples - 1;
               
                signalN = signal;

                isBefore = iStart <= 0;
                n_nanpad_presignal = 0;
                if any(isBefore)

                    n_nanpad_presignal = iStart(find(~isBefore,1,'first')) - iStart(1);
                    signalN = [nan(nChannels, n_nanpad_presignal), signalN];

                end

                isAfter = iEnd > nSample;
                if any(isAfter)

                    n_nanpad_postsignal = iEnd(end) - nSample;
                    signalN = [signalN, nan(nChannels, n_nanpad_postsignal)];

                end

                ep = nan([nTrials, nChannels, nEpochSamples]);

                for iEp = 1: nTrials

                    for iCh = 1:nChannels

                        ep(iEp, iCh, :) = signalN(iCh, n_nanpad_presignal + (iStart(iEp):iEnd(iEp)));

                    end

                end


            end

            
        end

        

    end
    methods

        function ep = retrieve(eTbl, varargin)

            ep = ns.RetrievedEpochs(eTbl, varargin{:});

        end


        

    end

    methods % GET Methods

        function c = get.C(eTbl)

            c = ns.C & eTbl;

        end
        function ch = get.channels(eTbl)

            ch = eTbl.C.channels;
            
        end
        
        function trials = get.trials(eTbl)

            eTbl.trials_.update(eTbl);

            trials = eTbl.trials_.value;


        end
    end

    methods (Static)

        function ep = detrend(ep, varargin)

            n_dim = ndims(ep);
            ep = permute(ep, circshift(1:n_dim,1));
            ep = detrend(ep,1,'omitnan', varargin{:});
            ep = permute(ep, circshift(1:n_dim,-1));

        end

        function ep = baseline(ep, isEpWin)

            base = mean(ep(:,:,isEpWin), ndims(ep), 'omitnan');
            ep = ep - base;

        end

        function ep = rereference(ep, chanLocs, varargin)

            varargin = [varargin, 'channel_locations', chanLocs];

            ep = ns.rereference(ep, varargin{:});

        end

    end

end

