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

    methods (Access=public)
        function T = find(tbl,value,pv)
            % Find components matching a query. See ns.Label.find for
            % examples
            % Delegates to ns.Label.findInTable so the logic is maintained in one place.
            arguments
                tbl (1,1) ns.LabelSession
                value  (1,:)  {mustBeA(value,["string" "double"])} % One or more values to look for
                pv.op (1,:)  {mustBeA(pv.op,["function_handle" "string" "char"])} = function_handle.empty  % Operator to use. Defaults to == for string and > for numeric                
                pv.findExtra (1,1) logical =false 
            end
            % Pass to static that is also used by LabelSession
            pv = namedargs2cell(pv);
            T = ns.Label.findInTable(fetchtable(tbl * ns.LabelParm, '*'), value,pv{:});
        end
    end

    methods (Access=protected)
        function makeTuples(tbl, key)
            parms = fetch1(ns.LabelParm & key, 'parms');

            % All experiments in this session for the relevant ctag
            sessionKey = struct('subject', key.subject, ...
                                'session_date', key.session_date, ...
                                'ctag',         key.ctag);
            cKeys = fetch(ns.C * proj(ns.Experiment) & sessionKey);
            assert(~isempty(cKeys), ...
                'No ns.C data found for ctag "%s" in session %s / %s.', ...
                key.ctag, key.subject, string(key.session_date));

            % Load all EEG datasets for this session, then delegate to the
            % shared computation in ns.Label.computeLabels.
            EEGs = cell(1, numel(cKeys));
            for i = 1:numel(cKeys)
                EEGs{i} = ephys.eeglab.dataset(cKeys(i), data=key.ctag, itag=key.itag);
            end
            [q, extra] = ns.Label.computeLabels(parms, EEGs, cKeys);
            insert(tbl, mergestruct(key, struct('q', q, 'extra', extra)));
        end
    end

end
