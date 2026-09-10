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
            % Sessions where ns.C is FULLY populated for the ctag specified
            % in a session-scoped IcaParm.  A session is only eligible once
            % every experiment that is in ns.C's keySource for that ctag
            % actually has a populated ns.C row.

            sessionParms = proj(ns.IcaParm & 'session=1', 'itag', 'ctag');
            % Candidate sessions: at least one ns.C row exists for the ctag.
            candidates = proj(ns.Session) * sessionParms ...
                & proj(ns.C, 'subject', 'session_date', 'ctag');

            if ~exists(candidates)
                v = candidates & 'FALSE';
                return;
            end

            % Single-query completeness check:
            % Restrict ns.C's keySource to the candidate (session, ctag) pairs,
            % then subtract what is already populated.  Any remaining rows identify
            % sessions that still have unpopulated entries.
            candSessionCtag = proj(candidates, 'subject', 'session_date', 'ctag');
            missing = (ns.C().keySource & candSessionCtag) - (ns.C & candSessionCtag);

            if exists(missing)
                fprintf(2,'ns.C incomplete for the following sessions (sessions skipped for ICA):\n')
                disp(ns.Session & missing )
                % Exclude incomplete sessions from the keySource.
                v = candidates - missing;
            else
                v = candidates;
            end
        end
    end

    methods (Access=public)

        function plot(tbl, pv)
            arguments
                tbl (1,1) ns.IcaSession
                pv.comp (1,:) double {mustBeInteger,mustBePositive} = 1:12
                pv.labels (1,:) string = ""  % Which ns.Label ltag to use for labeling components. Defaults to all
                pv.find = string.empty  % Optional string to filter components by their q in ns.Label
                pv.tilesPerFigure (1,1) double = 24
                pv.etaDisplay (1,1) string {mustBeMember(pv.etaDisplay, ["button", "inline", "none"])} = "button"
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

                % Fetch labels from ns.LabelSession (computed once for the whole session).
                labelRelvar = ns.LabelSession * ns.LabelParm & this;
                if pv.labels~=""
                    labelRelvar = labelRelvar & struct('ltag', cellstr(pv.labels)');
                end
                labels = fetch(labelRelvar, 'q', 'extra', 'parms');

                compsToPlot = pv.comp;
                if ~isempty(pv.find)
                    if ~isempty(labels)
                        foundComps = find(ns.LabelSession & this, pv.find);
                        if isempty(foundComps)
                            warning('No components found matching find instruction.');
                            continue;
                        else
                            compsToPlot = [foundComps.components{:}];
                        end
                    else
                        warning('Cannot filter by find: no ns.LabelSession rows found for this session ICA.');
                    end
                end

                expName = sprintf('%s @ %s | %s/%s', this.subject, this.session_date, this.ctag, this.itag);
                ns.Ica.plotComponents(this.winverse, this.variance, chanlocs, labels, compsToPlot, expName, pv.tilesPerFigure, pv.etaDisplay);
            end
        end

    end

    methods (Access=protected)
        function makeTuples(tbl, key)
            % Load all EEG datasets for this session and ctag, concatenate,
            % run ICA once, and store the weights.
            % DataJoint's jobs table ensures this runs exactly once per key.

            icaParms = fetch1(ns.IcaParm & key, 'parms');

            cKeys = fetch(ns.C & key & struct('ctag', key.ctag));


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