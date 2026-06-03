%{
# Labels for independent components 
-> ns.Ica
-> ns.LabelParm
---
q =NULL : longblob  # method-specific data used for labeling (e.g. ICLabel most likely class or ETA-based z-scores)
extra =NULL : longblob # method-specific extra data (e.g. ICLabel full class probabilities or ETA-based event-triggered averages)
%}
%
% After running ICA (see ns.Ica and ns.IcaParm), we want to label
% components for their subsequent use in epoching (ns.Epoch) or to remove
% components during preprocessing.
%
% This table can be filled by different methods labeling the components.
% For instance the eeglab iclabel function will fill q with the highest
% probablity label for each component, and extra with the full probability
% matrix. See ns.Ica/plot() and ns.LabelParm how to set this up for
% populate.
%
% The eta method identifies labels by determining an event triggered average (eta)
% for a (set of) specified events in a plugin and assessing whether that
% event led to a signficant change in the componennt.  This is, for
% instance, a good way to identify blink related components.
% This methods stores a z-score per component; higher z-scores suggest that
% the event really drove the component. (z-score is the difference between
% the component activity at the time of the event compared to a baseline
% period just before).The extra field stores the baseline corrected event
% triggered averages for each component.
%
% See also ns.LabelParm
%
% BK - May 2026.

classdef Label <  dj.Computed & dj.DJInstance

    properties (Dependent)
        keySource
    end

    methods
        function v = get.keySource(~)
            % Produce (Ica, LabelParm) pairs respecting paradigm restrictions.
            %
            % For each ltag:
            %   - No rows in LabelParmParadigm → applies to all paradigms
            %   - Has rows in LabelParmParadigm → applies only to those paradigms
            %
            % Per-experiment ICA (IcaParm.session=0): cross-product of Ica x
            % LabelParm filtered by paradigm, as before.
            %
            % Session ICA (IcaParm.session=1): gated on ns.LabelSession being
            % populated so that DataJoint does not attempt per-experiment labels
            % before the session-level labels exist.
            allParms = fetch(ns.LabelParm, 'ltag');

            if isempty(allParms)
                v = ns.Ica * proj(ns.Experiment, 'paradigm') * ns.LabelParm & 'FALSE';
                return;
            end

            parmClauses = cell(numel(allParms), 1);
            for p = 1:numel(allParms)
                ltag = allParms(p).ltag;
                paradigms = fetchn(ns.LabelParmParadigm & struct('ltag', ltag), 'name');
                if isempty(paradigms)
                    % Unrestricted: applies to every paradigm
                    parmClauses{p} = sprintf('(ltag="%s")', ltag);
                else
                    % Restricted: only for the listed paradigms
                    paradigmClauses = cellfun(@(pg) sprintf('paradigm="%s"', pg), ...
                        paradigms, 'UniformOutput', false);
                    parmClauses{p} = sprintf('(ltag="%s" AND (%s))', ltag, ...
                        strjoin(paradigmClauses, ' OR '));
                end
            end

            combinedWhere = ['(' strjoin(parmClauses, ' OR ') ')'];

            % Non-session ICA: paradigm-filtered.
            nonSessionPart = (ns.Ica * proj(ns.Experiment, 'paradigm') * ...
                proj(ns.IcaParm & 'session=0') * ...
                ns.LabelParm) & combinedWhere;

            % Session ICA: one Label row per experiment, driven by LabelSession.
            % The natural join on (subject, session_date, itag, ctag, ltag)
            % produces one row per Ica entry once the session label exists.
            sessionPart = ns.Ica * ...
                proj(ns.IcaParm & 'session=1') * ...
                proj(ns.LabelSession, 'ltag');

            % DataJoint MATLAB has no union (+) operator.  Instead, fetch
            % both key sets and combine as a struct-array restriction on the
            % base relation.
            nonSessionKeys = fetch(nonSessionPart);   % returns Label PK fields only
            sessionKeys    = fetch(sessionPart);       % returns Label PK fields only
            allKeys = [nonSessionKeys; sessionKeys];
            if isempty(allKeys)
                v = ns.Ica * ns.LabelParm & 'FALSE';
            else
                v = ns.Ica * ns.LabelParm & allKeys;
            end
        end
    end


    methods (Access=public)
        function T = find(tbl, q)
            % Find components matching a query string (for iclabel) or numeric threshold (for eta method).
            arguments
                tbl (1,1) ns.Label
                q (1,:)
            end
            T = ns.Label.findInTable(fetchtable(tbl * ns.LabelParm, '*'), q);
        end
    end

    methods (Static)
        function T = findInTable(T, q)
            % Core find logic shared by ns.Label.find and ns.LabelSession.find.
            % T must be a table with columns 'parms' and 'q' (from a *ns.LabelParm join).
            T = addvars(T, cell(height(T), 1), 'NewVariableNames', 'components');
            out = false(height(T),1);
            for tpl = 1:height(T)
                if isstruct(T.parms)
                    method = T.parms(tpl).method;
                else
                    method = T.parms{tpl}.method;
                end
                id = sprintf("%s@%sT%s",T.subject(tpl),T.session_date(tpl),T.starttime(tpl));
                switch upper(method)
                    case 'ICLABEL'
                        if ~iscellstr(q) && ~isstring(q) && ~ischar(q)
                            out(tpl) =true;
                            continue;
                        end
                        comp = find(contains(T.q{tpl}, q, 'IgnoreCase', true));
                        if isempty(comp)
                            fprintf('No components matching "%s" in %s.\n', q,id);
                        else
                            fprintf('Components matching "%s" in %s: %s\n', q, id,strjoin(string(comp), ', '));
                        end
                    case {'ETA','SPEARMAN'}
                        if ~isnumeric(q)
                            out(tpl) =true;
                            continue;
                        end
                        comp = find(abs(T.q{tpl}) > q);
                        if isempty(comp)
                            fprintf('No components with |q|>=%.2f in %s.\n', q,id);
                        end
                    otherwise
                        error('No find implemented for labeling method: %s', method);
                end
                T{tpl, 'components'} = {comp};
            end
            T(out,:) = [];
        end
    end




    methods (Static)
        function [q, extra] = computeLabels(parms, EEGs, cKeys)
            % Compute ICA component labels for one or more EEG datasets.
            %
            % parms   - struct from ns.LabelParm.parms
            % EEGs    - 1xN cell array of EEG structs (ICA loaded, icaact populated)
            % cKeys   - N-element struct array of C keys, one per EEG.
            %           Used by SPEARMAN to fetch the auxiliary channel signal.
            %           Pass [] for ICLABEL and ETA.
            %
            % Returns q and extra with the same semantics as ns.Label.q / extra.

            switch upper(parms.method)

                case 'ICLABEL'
                    % ICLabel depends only on component topographies (winverse),
                    % which are session-invariant, so one EEG suffices.
                    assert(exist('iclabel', 'file') == 2, ...
                        'EEGLAB with iclabel function is required for method ''iclabel''.');
                    EEG   = iclabel(EEGs{1});
                    extra = EEG.etc.ic_classification.ICLabel.classifications;
                    [~, ix] = max(extra, [], 2);
                    q = {EEG.etc.ic_classification.ICLabel.classes(ix)};

                case 'SPEARMAN'
                    allAct    = [];
                    allSignal = [];
                    for i = 1:numel(EEGs)
                        EEG = EEGs{i};
                        % Strip fields that are not part of the C primary key
                        % so the restriction works regardless of caller context.
                        fieldsToStrip = intersect(fieldnames(cKeys(i)), {'ctag','filename','itag','ltag'});
                        expKey = rmfield(cKeys(i), fieldsToStrip);
                        c = ns.C & expKey & struct('ctag', parms.ctag);
                        if ~exists(c)
                            fprintf(2,"No %s data in %s@%sT%s\n",parms.ctag,expKey.subject,expKey.session_date,expKey.starttime)
                        else
                            tplC = fetch(c * ns.CChannel & struct('name', parms.channel), '*');
                            eegSampleTime = polyval(EEG.etc.neurostim.clockParms, EEG.times/1000)';
                            cSampleTime   = linspace(tplC.time(1), tplC.time(2), tplC.time(3))';
                            cSignal       = interp1(cSampleTime, tplC.signal, eegSampleTime);
                            cSignal(isnan(cSignal)) = tplC.min;
                            allAct    = [allAct,    EEG.icaact]; %#ok<AGROW>
                            allSignal = [allSignal; cSignal];    %#ok<AGROW>
                        end
                    end
                    if isempty(allAct)
                        q= [];
                        extra = [];
                    else
                        [q, extra] = corr(allAct', allSignal, 'Type', 'Spearman');
                    end

                case 'ETA'
                    parms.events = string(parms.events);
                    EEG0     = EEGs{1};
                    windowMs = double(parms.window(:)');
                    baseMs   = double(parms.baseline(:)');
                    sampleOffsets = round(windowMs(1) * EEG0.srate / 1000):round(windowMs(2) * EEG0.srate / 1000);
                    nSamps   = numel(sampleOffsets);
                    timesMs  = sampleOffsets / EEG0.srate * 1000;
                    nComps   = size(EEG0.icaact, 1);
                    allEpochs = zeros(nComps, nSamps, 0);

                    for i = 1:numel(EEGs)
                        EEG = ephys.eeglab.addEvents(EEGs{i}, parms.plugin, parms.events);
                        allEventTypes = erase(lower(string({EEG.event.type})), "'");
                        keep      = ismember(allEventTypes, parms.events);
                        latencies = round([EEG.event(keep).latency]);
                        ok        = (latencies + sampleOffsets(1)) >= 1 & (latencies + sampleOffsets(end)) <= EEG.pnts;
                        latencies = latencies(ok);
                        if isempty(latencies), continue; end
                        idx    = (sampleOffsets + latencies(:))';
                        epochs = reshape(EEG.icaact(:, idx), [nComps nSamps numel(latencies)]);
                        allEpochs = cat(3, allEpochs, epochs);
                    end
                    nEvents = size(allEpochs,3);
                    if isfield(parms,'minNrEvents')
                        minNrEvents = parms.minNrEvents;
                    else
                        minNrEvents = 10;
                    end
                    if nEvents < minNrEvents
                        warning('ns:Label:nrEvents', 'Only %d %s events found. Skipping.', nEvents,strjoin(parms.events));
                        q = []; extra = [];
                    else
                        baseWin = timesMs >= baseMs(1) & timesMs <= baseMs(2);
                        % Determine whether there are timepoints that
                        % differ significantly from the baseline by doing
                        % consecutive t-test, followed by FDR correction for the number of components
                        % and the number of timepoints. 
                        % The q value is -log10(pCorrected). 1 means pCorrected = 0.05.
                        % The extra is the event triggered average.

                        % Baseline correction
                        mu_base = mean(allEpochs(:, baseWin,  :), 2);
                        allEpochs = (allEpochs - mu_base);
                        extra=  mean(allEpochs, 3, 'omitmissing');

                        % Test each baseline corrected component at each sample against 0 
                        allEpochs = allEpochs(:,~baseWin,:);
                        nSamps =sum(~baseWin);
                        p = nan(nComps,nSamps);
                        t = nan(nComps,nSamps);
                        for c = 1:size(allEpochs,1)                            
                                [~, p(c,:), ~, stats] = ttest(squeeze(allEpochs(c,:,:))', 0);
                                t(c, :) = [stats.tstat];                            
                        end
                        % FDR correction at alpha=0.05
                        alpha =1/(nComps+1);
                        [pFdr] = fdr(p(:), alpha);                        
                        [pMin,timePoint] = min(p,[],2);
                        samples = reshape(allEpochs,[nComps*nSamps nEvents]);
                        ix = sub2ind([nComps nSamps],(1:nComps)',timePoint);
                        samples = reshape(samples(ix,:),[nComps nEvents]);
                        q = abs(mean(samples,2,"omitmissing"))./std(samples,0,2,"omitmissing");
                        q(pMin>pFdr) =0;
                        
                    end

                otherwise
                    error('Unknown labeling method: %s', parms.method);
            end
        end
    end

    methods (Access=protected)
        function makeTuples(tbl, key)
            % For session-scoped ICA, labels are computed once in ns.LabelSession;
            % copy them here rather than recomputing per experiment.
            isSession = fetch1(ns.IcaParm & key, 'session');
            if isSession
                sessionKey = struct('subject', key.subject, 'session_date', key.session_date, ...
                    'itag', key.itag, 'ctag', key.ctag, 'ltag', key.ltag);
                src = fetch1(ns.LabelSession & sessionKey, 'q', 'extra');
                insert(tbl, mergestruct(key, src));
                return;
            end

            parms = fetch1(ns.LabelParm & key, 'parms');
            EEG   = ephys.eeglab.dataset(key, data=key.ctag, itag=key.itag);
            [q, extra] = ns.Label.computeLabels(parms, {EEG}, key);
            insert(tbl, mergestruct(key, struct('q', q, 'extra', extra)));
        end
    end


end



