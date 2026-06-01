%{
# Labels for ICA components computed at the session level
-> ns.IcaSession
-> ns.LabelParm
---
q     =NULL : longblob   # method-specific data (same semantics as ns.Label.q)
extra =NULL : longblob   # method-specific extra data (same semantics as ns.Label.extra)
%}
%
% This table stores component labels computed once per session for a
% session-scoped ICA decomposition (ns.IcaSession).  It mirrors ns.Label
% but operates on the full merged session dataset rather than on individual
% experiments, so that labels such as blink/saccade ETAs or Spearman
% correlations benefit from the larger sample size.
%
% For ICLABEL the labels depend only on component topographies (winverse),
% which are identical across all experiments in the session, so only one
% experiment needs to be loaded.  For ETA and SPEARMAN the ICA activations
% and the auxiliary signal are concatenated across all experiments in the
% session before computing the statistics.
%
% Populate this table BEFORE ns.Label when using session-scoped ICA.
% ns.Label.makeTuples detects IcaParm.session=1 and copies q/extra from
% this table rather than recomputing per experiment.
%
% See also ns.Label, ns.LabelParm, ns.IcaSession
%
% BK - May 2026

classdef LabelSession < dj.Computed & dj.DJInstance

    properties (Dependent)
        keySource
    end

    methods
        function v = get.keySource(~)
            v = ns.IcaSession * ns.LabelParm;
        end
    end

    methods (Access=protected)
        function makeTuples(tbl, key)
            parms = fetch1(ns.LabelParm & key, 'parms');

            % All experiments in this session for the relevant ctag
            sessionKey = struct('subject',      key.subject, ...
                                'session_date', key.session_date, ...
                                'ctag',         key.ctag);
            cKeys = fetch(ns.C * proj(ns.Experiment) & sessionKey);
            assert(~isempty(cKeys), ...
                'No ns.C data found for ctag "%s" in session %s / %s.', ...
                key.ctag, key.subject, string(key.session_date));

            switch upper(parms.method)

                case 'ICLABEL'
                    % ICLabel depends only on topographies, which are
                    % session-invariant.  Load one experiment.
                    assert(exist('iclabel', 'file') == 2, ...
                        'EEGLAB with iclabel is required for method ''iclabel''.');
                    EEG   = ephys.eeglab.dataset(cKeys(1), data=key.ctag, itag=key.itag);
                    EEG   = iclabel(EEG);
                    extra = EEG.etc.ic_classification.ICLabel.classifications;
                    [~, ix] = max(extra, [], 2);
                    q = {EEG.etc.ic_classification.ICLabel.classes(ix)};

                case 'SPEARMAN'
                    % Concatenate ICA activations and the auxiliary signal
                    % across all experiments in the session.
                    allAct    = [];
                    allSignal = [];
                    for i = 1:numel(cKeys)
                        EEG = ephys.eeglab.dataset(cKeys(i), data=key.ctag, itag=key.itag);
                        c    = ns.C & rmfield(cKeys(i), ["ctag" "filename"]) & struct('ctag', parms.ctag);
                        tplC = fetch(c * ns.CChannel & struct('name', parms.channel), '*');
                        eegSampleTime = polyval(EEG.etc.neurostim.clockParms, EEG.times/1000)';
                        cSampleTime   = linspace(tplC.time(1), tplC.time(2), tplC.time(3))';
                        cSignal       = interp1(cSampleTime, tplC.signal, eegSampleTime);
                        cSignal(isnan(cSignal)) = tplC.min;
                        allAct    = [allAct,    EEG.icaact]; %#ok<AGROW>
                        allSignal = [allSignal; cSignal];    %#ok<AGROW>
                    end
                    [q, extra] = corr(allAct', allSignal, 'Type', 'Spearman');

                case 'ETA'
                    parms.events = string(parms.events);
                    EEG0     = ephys.eeglab.dataset(cKeys(1), data=key.ctag, itag=key.itag);
                    preSamp  = round(parms.window * EEG0.srate / 1000);
                    postSamp = preSamp;
                    nSamps   = preSamp + postSamp + 1;
                    timesMs  = (-preSamp:postSamp) / EEG0.srate * 1000;
                    nComps   = size(EEG0.icaact, 1);
                    allEpochs = zeros(nComps, nSamps, 0);

                    for i = 1:numel(cKeys)
                        EEG = ephys.eeglab.dataset(cKeys(i), data=key.ctag, itag=key.itag);
                        EEG = ephys.eeglab.addEvents(EEG, parms.plugin, parms.events);
                        allEventTypes = erase(lower(string({EEG.event.type})), "'");
                        keep      = ismember(allEventTypes, parms.events);
                        latencies = round([EEG.event(keep).latency]);
                        ok        = latencies > preSamp & latencies <= (EEG.pnts - postSamp);
                        latencies = latencies(ok);
                        if isempty(latencies), continue; end
                        idx    = ((-preSamp:postSamp) + latencies(:))';
                        epochs = reshape(EEG.icaact(:, idx), [nComps nSamps numel(latencies)]);
                        allEpochs = cat(3, allEpochs, epochs);
                    end

                    if size(allEpochs, 3) == 0
                        warning('ns:LabelSession:noEvents', ...
                            'No %s events found across session %s / %s.', ...
                            strjoin(parms.events), key.subject, string(key.session_date));
                        q = []; extra = [];
                    else
                        baseWin = timesMs >= -parms.window & timesMs < -parms.window + parms.baseline;
                        mu_base = mean(allEpochs(:, baseWin,  :), [2 3]);
                        sd_base = std( allEpochs(:, baseWin,  :), 0, [2 3]);
                        mu_resp = mean(allEpochs(:, ~baseWin, :), [2 3]);
                        q       = (mu_resp - mu_base) ./ sd_base;
                        extra   = mean(allEpochs, 3, 'omitmissing') - mu_base;
                    end

                otherwise
                    error('Unknown labeling method: %s', parms.method);
            end

            insert(tbl, mergestruct(key, struct('q', q, 'extra', extra)));
        end
    end

end
