%{
# Parameters that define a set of conditions (a Dimension) for a paradigm.
dimension : varchar(32) # Name of this dimension
paradigm : varchar(255) # Paradigm for which this dimension is defined (matches ns.Experiment.paradigm)
---
parms  : longblob  # Parameters that define the condition
description = NULL : varchar(1024) # Description
%}
% Each row defines how a named dimension is computed for a specific paradigm.
% The compound primary key (dimension, paradigm) allows the same conceptual
% dimension to have different plugin/parameter definitions in different paradigms.
%
% EXAMPLE
%  For the 'SSVEP' paradigm, define a dimension called 'frequency' driven by
%  the 'frequency' property of the 'flicker' plugin:
%
%   tpl = struct('paradigm', 'SSVEP', ...
%               'dimension', 'frequency', ...
%               'parms', struct('plg', 'flicker', 'prm', 'frequency'));
%   insert(ns.DimensionParm, tpl);
%
% EXAMPLE 2
%  Optional parms fields restrict or exclude trials: 
% (Note the double {{}}; those are needed to get a scalar parms struct
%   parms = struct('plg', 'flicker', 'prm', 'frequency', ...
%                  'restrict', {{'cic','block',2}}, ...
%                  'exclude',  {{'flicker','duration',0}});
%
% See also ns.Dimension

classdef DimensionParm < dj.Lookup


    methods (Access=public)
        function insert(self, tuples, varargin)
            % Overload insert to validate parms and skip/error on duplicates.
            for t = 1:numel(tuples)
               assert(isscalar(string(tuples(t).paradigm)),"DimensionParms can list a single paradigm only.");
               assert(isscalar(tuples(t).parms),"Parms must be a scalar struct. Use double curly brackets for cell arrays (restrict and exclude) {{plg,prm,val}}");
                pv = namedargs2cell(tuples(t).parms);
                tuples(t).parms = ns.DimensionParm.validateParms( ...
                    string(tuples(t).paradigm), pv{:});
            end
            tuples = makeMymSafe(tuples);

            % Check for existing (dimension, paradigm) pairs
            isExisting = false(numel(tuples), 1);
            for t = 1:numel(tuples)
                key = struct('dimension', tuples(t).dimension, ...
                             'paradigm',  tuples(t).paradigm);
                existing = fetch(self & key, 'parms');
                if ~isempty(existing)
                    if isequal(existing.parms, tuples(t).parms)
                        fprintf('DimensionParm (''%s'', ''%s'') already exists with identical parms; skipping.\n', ...
                            tuples(t).dimension, tuples(t).paradigm);
                        isExisting(t) = true;
                    else
                        error('ns:DimensionParm:conflictingParms', ...
                            'DimensionParm (''%s'', ''%s'') exists with different parms. Cannot insert.', ...
                            tuples(t).dimension, tuples(t).paradigm);
                    end
                end
            end
            if any(~isExisting)
                insert@dj.Lookup(self, tuples(~isExisting), varargin{:})
            end
        end

        function validate(tbl, opts)
            arguments
                tbl (1,1) ns.DimensionParm
                opts.full (1,1) logical = false
            end
            % Validate parms for each (dimension, paradigm) row.
            for tpl = fetch(tbl, '*')'
                pv = namedargs2cell(tpl.parms);
                optsPv = namedargs2cell(opts);
                ns.DimensionParm.validateParms(string(tpl.paradigm), pv{:}, optsPv{:});
            end
        end
    end


    methods (Static)
        function pv = validateParms(paradigm, pv, opts)
            arguments
                paradigm (1,:) string = ""  % paradigm name for validation context
                pv.plg (1,:) string  =""
                pv.prm (1,:) string =""
                pv.fun (1,1) string = ""  % Name of a function that converts plg.prm into a number that defines the condition (default is to do nothing; the value determines the condition)
                pv.restrict (1,:) cell = {} % Define only in this subset of trials by specifying a set of allowed values for a plugin parameter
                pv.exclude (1,:) cell = {} % Define only in this subset of trials by specifying a set of disallowed values for a plugin parameter
                pv.left (1,1) double = 0 % Reduce the names to this number of chars from the left. 0 means the full length
                pv.nameValueOnly (1,1) = false  % Set to true to define condition names based o the prm values alone (and not their name).
                pv.atTrialTime (1,1) = 0 % By default a dimension is defined by the parameter value at the start of the trial. Set this to Inf to use the value at the end of the trial (or any other trial time).
                opts.full (1,1) logical =false;
            end

            assert(numel(pv.plg)==numel(pv.prm),'The dimension parm must define a parameter (.prm) for each plugin (.plg)');

            for thisParadigm = paradigm
                expts = proj(ns.Experiment & struct('paradigm',thisParadigm));
                if count(expts)==0
                    fprintf('No experiments with %s paradigm. Cannot validate.\n',thisParadigm);
                    continue;
                end
                % Check that plugins exist and have the parameter
                nrPlgsAndPrms = numel(pv.plg);
                for p = 1:nrPlgsAndPrms
                    thisPlg = ns.Plugin & struct('plugin_name',pv.plg(p)) & expts;
                    nrPlgs = count(thisPlg);
                    if nrPlgs==0
                        fprintf('No plugins named %s. Could be ok if you have not yet added any experiments with paradigm %s.\n',pv.plg(p),thisParadigm);
                        continue;
                    else
                        % Check that each plugin has the prm
                        thisPlgWithPrm= ns.PluginParameter & thisPlg & struct('property_name',pv.prm(p));
                        nrWithPrm = count(thisPlgWithPrm);
                        if (nrWithPrm~=nrPlgs)
                            fprintf('%d plugins named %s, but only %d have the property named %s\n. These do not:',nrPlgs,pv.plg(p),nrWithPrm,pv.prm(p));
                        end
                        if opts.full
                            % Could add checking of the actual values.
                        end
                    end

                end
            end
        end



    end
end