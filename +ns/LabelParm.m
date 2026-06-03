%{
# Parameter sets for determining labels for ICA 
ltag         :  varchar(32)      # A unique name for this label definition
---
description  :  varchar(1024)    # Short description
parms        :  longblob         # struct with Label-specific options
%}
%
% EXAMPLE: Apply to all paradigms (warns)
%   tpl = struct('ltag','iclabel', ...
%                'description','Label ICA components using EEGLAB ICLabel', ...
%                'parms',struct('method','iclabel','version','default'));
%   insert(ns.LabelParm, tpl)
%
% EXAMPLE: Restrict to specific paradigms
%   insert(ns.LabelParm, tpl, paradigm=["bandit" "ssvep"])
%
% EXAMPLE: Add paradigms to an existing LabelParm
%   ns.LabelParm.addToParadigm('eye', "ssvep")
%
% EXAMPLE: tpl = struct('ltag','blink', ...
%               'description','Label ICA components using Eyelink blink events', ...
%               'parms',struct('method','eta', ...
%                                'window',200,... % Extract a window [-200 +200] around each event
%                                'baseline',100,...  % Define the first 100 ms of the window as baseline
%                                'plugin','edf', ...
%                                'events','startblink'));
% See also ns.Label, ns.LabelParmParadigm
%
% BK - May 2026
classdef LabelParm < dj.Lookup

    methods
        function insert(self, tuples, pv)
            % Validate then insert into LabelParm.
            % Optionally restrict to specific paradigms.
            %
            % USAGE
            %   insert(ns.LabelParm, tpl)
            %       Inserts tpl and warns that it will apply to all paradigms.
            %   insert(ns.LabelParm, tpl, paradigm=["a" "b"])
            %       Inserts tpl and adds rows to LabelParmParadigm for each paradigm.
            arguments
                self   (1,1) ns.LabelParm
                tuples (1,:) struct
                pv.paradigm (1,:) string = string.empty
            end
            ns.LabelParm.validate(tuples);            
            insert@dj.Lookup(self, makeMymSafe(tuples));
            if isempty(pv.paradigm)
                warning('ns:LabelParm:allParadigms', ...
                    'No paradigm specified for ltag "%s": it will apply to ALL paradigms.', ...
                    tuples.ltag);
            else
                ns.LabelParm.addToParadigm(tuples.ltag, pv.paradigm);
            end
        end
    end

    methods (Static)
        function addToParadigm(ltag, paradigmName)
            % Add an existing LabelParm entry to one or more paradigms.
            % paradigmName entries containing '%' are treated as SQL LIKE
            % patterns and expanded to all matching rows in ns.Paradigm.
            %
            % USAGE
            %   ns.LabelParm.addToParadigm('blink', 'myParadigm')
            %   ns.LabelParm.addToParadigm('blink', ["paradigm1" "paradigm2"])
            %   ns.LabelParm.addToParadigm('blink', "kEEG%")   % wildcard
            arguments
                ltag         (1,1) string
                paradigmName (1,:) string
            end
            resolved = string.empty(1,0);
            for i = 1:numel(paradigmName)
                if contains(paradigmName(i), '%')
                    matches = fetchn(ns.Paradigm & sprintf('name LIKE "%s"', paradigmName(i)), 'name');
                    assert(~isempty(matches), ...
                        'No paradigms in ns.Paradigm match pattern "%s".', paradigmName(i));
                    resolved = [resolved, string(matches(:)')]; %#ok<AGROW>
                else
                    resolved(end+1) = paradigmName(i); %#ok<AGROW>
                end
            end
            for i = 1:numel(resolved)
                tpl = struct('ltag', char(ltag), 'name', char(resolved(i)));
                insert(ns.LabelParmParadigm, tpl);
            end
        end
    end

    methods (Static, Access = protected)
        function validate(tuples)
            % Validate LabelParm tuples before insertion.
            % Add field checks and constraints here.
            for tpl = tuples
                switch tpl.parms.method
                    case "iclabel"
                        if ~exist('iclabel', 'file')
                            warning('ICLabel method specified but iclabel.m not found on path.');
                        end
                        assert(isfield(tpl.parms, 'version'), 'Missing "version" field in iclabel parm.');
                    case "eta"
                        assert(isfield(tpl.parms, 'plugin'), ...
                            'eta method requires a "plugin" field (name of the Neurostim plugin containing events).');
                        assert(isfield(tpl.parms, 'events'), ...
                            'eta method requires an "events" field (string or string array of event names).');
                        assert(isfield(tpl.parms, 'window'), ...
                            'eta method requires a "window" field ([start stop] in ms relative to the event).');                        
                        assert(isfield(tpl.parms, 'baseline'), ...
                            'eta method requires a "baseline" field ([start stop] in ms relative to the event).');
                        assert(isnumeric(tpl.parms.window) && numel(tpl.parms.window)==2 , ...
                            '"window" must be  [ start stop] (ms).');
                        assert(isnumeric(tpl.parms.baseline) && numel(tpl.parms.baseline)==2 , ...
                            '"baseline" must be [start stop] (ms).');
                        assert(tpl.parms.window(1) <= tpl.parms.window(2), ...
                            '"window" must satisfy start <= stop.');
                        assert(tpl.parms.baseline(1) <= tpl.parms.baseline(2), ...
                            '"baseline" must satisfy start <= stop.');
                        assert(tpl.parms.baseline(1) >= tpl.parms.window(1) && tpl.parms.baseline(2) <= tpl.parms.window(2), ...
                            '"baseline" must lie within "window".');
                    case "spearman"
                        assert(isfield(tpl.parms,"ctag"),...
                            'spearman method requires a "ctag" field specifying which ns.C data to use for regression.');
                        assert(isfield(tpl.parms,"channel"),...      
                            'spearman method requires a "channel" field specifying which channel in ns.C to use for regression.');                    

                    otherwise
                        error('Unsupported labeling method: %s', tpl.parms.method);
                end
            end
        end
    end

end
