%{
# ICA computed across all experiments in a session (one decomposition per session/itag)
-> ns.Session
-> ns.IcaParm
---
nrcomponents : int       # Number of components
nrsamples    : int       # Number of samples used for ICA
winverse     : longblob  # Inverse ICA weights (mixing matrix)
sphere       : longblob  # Sphering matrix
weights      : longblob  # ICA weights (unmixing)
chanlabels   : longblob  # Channel labels (strings) used for ICA - intersection across all experiments
variance     : blob      # Variance explained per component
%}
%
% This table stores a single ICA decomposition computed from the
% concatenation of all ns.C datasets in a session for a given ctag/itag.
% It is populated before ns.Ica when IcaParm.session=true.
% ns.Ica.makeTuples then copies the weights from this table into one
% ns.Ica row per experiment, so all downstream tables (ns.Label, ns.Epoch)
% remain unchanged.
%
% See also ns.Ica, ns.IcaParm
%
% BK - May 2026
classdef IcaSession < dj.Computed & dj.DJInstance

    properties (Dependent)
        keySource
    end

    methods
        function v = get.keySource(~)
            % All sessions that have at least one ns.C entry for the ctag
            % specified in a session-scoped IcaParm.
            sessionParms = proj(ns.IcaParm & 'session=1', 'itag', 'ctag');
            % Join Session with Experiment->C to verify data exists for ctag
            v = proj(ns.Session) * sessionParms ...
                & proj(ns.C * proj(ns.Experiment), 'subject', 'session_date', 'ctag');
        end
    end

    methods (Access=public)

        function plot(tbl, pv)
            arguments
                tbl (1,1) ns.IcaSession
                pv.comp (1,:) double {mustBeInteger,mustBePositive} = 1:12
                pv.labels (1,:) string = "iclabel"  % Which ns.Label ltag to use for labeling components
                pv.find = string.empty  % Optional string to filter components by their q in ns.Label
                pv.tilesPerFigure (1,1) double = 24
            end
            tpl = fetch(tbl, '*');
            assert(~isempty(tpl), 'No rows in ns.IcaSession to plot.');

            for i = 1:numel(tpl)
                this = tpl(i);
                sessionKey = struct('subject', this.subject, 'session_date', this.session_date);

                % Resolve channel locations: find any C experiment in session, filter to chanlabels.
                cKey = fetch(ns.C * proj(ns.Experiment) & sessionKey & struct('ctag', this.ctag), 'LIMIT 1');
                assert(~isempty(cKey), ...
                    'No ns.C data found for session %s/%s ctag=%s.', ...
                    this.subject, this.session_date, this.ctag);
                allChanlocs = [fetch(ns.CChannel & cKey(1), 'channelinfo').channelinfo];
                [~, idx] = ismember(this.chanlabels, {allChanlocs.labels});
                chanlocs = allChanlocs(idx(idx > 0));

                % Use any ns.Ica row for this session/itag to access labels
                % (all per-experiment rows share the same decomposition).
                icaRow = fetch(ns.Ica * proj(ns.Experiment) & sessionKey & struct('itag', this.itag), 'LIMIT 1');
                if ~isempty(icaRow)
                    labels = fetch(ns.Label * ns.LabelParm & icaRow(1) & struct('ltag', cellstr(pv.labels)'), 'q', 'extra', 'parms');
                else
                    labels = struct.empty;
                end

                compsToPlot = pv.comp;
                if ~isempty(pv.find)
                    if ~isempty(icaRow)
                        foundComps = find(ns.Label & icaRow(1), pv.find);
                        if isempty(foundComps)
                            warning('No components found matching find instruction. Plotting all components instead.');
                        else
                            compsToPlot = [foundComps.components{:}];
                        end
                    else
                        warning('Cannot filter by find: no ns.Ica/Label rows found for this session ICA.');
                    end
                end

                expName = sprintf('%s @ %s | %s/%s', this.subject, this.session_date, this.ctag, this.itag);
                ns.Ica.plotComponents(this.winverse, this.variance, chanlocs, labels, compsToPlot, expName, pv.tilesPerFigure);
            end
        end

    end

    methods (Access=protected)
        function makeTuples(tbl, key)
            % Load all EEG datasets for this session and ctag, concatenate,
            % run ICA once, and store the weights.
            % DataJoint's jobs table ensures this runs exactly once per key.

            icaParms = fetch1(ns.IcaParm & key, 'parms');

            % Find all experiments in this session that have C data for ctag
            sessionKey = rmfield(key, 'itag');
            cKeys = fetch(ns.C * proj(ns.Experiment) & sessionKey & struct('ctag', key.ctag));
            assert(~isempty(cKeys), ...
                'No ns.C data found for ctag "%s" in session %s/%s.', ...
                key.ctag, key.subject, key.session_date);

            % Load all datasets and find the intersection of surviving channel labels.
            % Experiments may have different bad channels removed during preprocessing,
            % so we can only run ICA on channels present in every experiment.
            EEGs = cell(1, numel(cKeys));
            for i = 1:numel(cKeys)
                EEGs{i} = ephys.eeglab.dataset(cKeys(i), data=key.ctag, itag="");
            end
            commonLabels = {EEGs{1}.chanlocs.labels};
            for i = 2:numel(EEGs)
                commonLabels = intersect(commonLabels, {EEGs{i}.chanlocs.labels}, 'stable');
            end
            assert(~isempty(commonLabels), ...
                'No common channels found across experiments in session %s/%s for ctag "%s".', ...
                key.subject, key.session_date, key.ctag);
            fprintf('Session ICA: using %d channels common to all %d experiments.\n', ...
                numel(commonLabels), numel(cKeys));

            % Reduce each dataset to the common channel set then merge
            EEG = pop_select(EEGs{1}, 'channel', commonLabels);
            for i = 2:numel(EEGs)
                EEGi = pop_select(EEGs{i}, 'channel', commonLabels);
                EEG  = pop_mergeset(EEG, EEGi);
            end

            % Optional high-pass filter before ICA
            if isfield(icaParms, 'filt')
                EEG = pop_eegfiltnew(EEG, 'hicutoff', icaParms.filt.hicutoff, ...
                                          'locutoff',  icaParms.filt.locutoff);
                icaParms = rmfield(icaParms, 'filt');
            end

            parms.eeglab.ica = icaParms;
            EEG = ephys.eeglab.preprocess(EEG, parms);

            varExplained = ephys.eeglab.icaVarianceExplained(EEG);

            tpl = key;
            tpl.nrcomponents = size(EEG.icaweights, 1);
            tpl.nrsamples    = EEG.pnts;
            tpl.winverse     = EEG.icawinv;
            tpl.sphere       = EEG.icasphere;
            tpl.weights      = EEG.icaweights;
            tpl.chanlabels   = commonLabels;  % channel labels used; remap to per-exp indices in getWeights
            tpl.variance     = varExplained;
            insert(tbl, tpl);
        end
    end

end