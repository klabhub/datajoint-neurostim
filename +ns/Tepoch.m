%{
# Transformed Epoch - a computation applied to a (set of) epochs.
-> ns.Epoch                     # The epochs that were transformed
-> ns.TepochParm                # Parameters used for the transformation
dependent       : varchar(64)   # The name of the dependent variable(s)
---
x : blob                        # The values of the independent variable
independent     : varchar(64)   #  The name of the independent variable(s), if multiple, concatenated with ':' in order
groups = NULL   : blob          # contains mappings for groups
%}
classdef Tepoch < dj.Computed & dj.DJInstance

    properties (SetAccess = protected)
        keySource = ns.TepochParm * (ns.Epoch & ns.EpochChannel)
    end
    methods (Access = protected)
        function makeTuples(tbl, key)
            % Apply a computation/transform to a collection of Epochs and
            % store as Tepoch.
            parms = fetch(ns.TepochParm & key,'*');

            % --- Select trials, conditions, and/or channels ---
            % TBD
            % --- Average by groups ---      
            eTbl = ns.Epoch & key;
            if isnan(parms.average)
                ecTbl = ns.EpochChannel & eTbl;
            else
                if ischar(parms.average), parms.average = {parms.average}; end
                [ecTbl, avg_groups] = find_groups_to_average_(tbl, key, parms.average{:});
            end
            % Restrict the epoch channels and trials if requested in the
            % parms
            if isempty(parms.window)
                % Use the same window as the epoch
                parms.window = fetch(ns.EpochParm & key,'window');
            end

            % Export groups
            gru_id = avg_groups.id;
            gru = rmfield(avg_groups, {'id','name'}); % must leave only the grouping variables
            gru_by = fieldnames(gru);
            compute_args = namedargs2cell(gru);
            n_gru = numel(gru_id);
            T = [];
            for iGru = 1:n_gru

                compute_argsN = compute_args;
                compute_argsN{2:2:end} = compute_argsN{2:2:end}{iGru};                

                % Compute - uses the ns.cache/compute function  (abstract
                % superclass).
                [tblN,dv,idv] = compute(ecTbl,parms.fun, compute_argsN{:}, average=gru_by, timeWindow=parms.window);
                
                % insert back the group variables (1st trial/channel)
                gru_col =  compute_argsN;
                gru_col{2:2:end} = gru_col{2:2:end}(1);
                gru_col = repelem(struct2table(struct(gru_col{:})), height(tblN),1);
                tblN.group = repelem(gru_id(iGru),height(tblN),1);
                tblN = [tblN, gru_col];
                T = [T; tblN];

                % I doubt if this work for a general case

            end           
           
            % Insert in the table
            %%
            x = table2cell(T(1,idv));
            x = cat(2,x{:});
            tpl = dj.struct.join(struct(independent = strjoin(idv,':'), x = x, dependent = cellstr(dv(:)), groups=avg_groups),key);

            insert(tbl,tpl);            

            dat_tbl = T(:,["channel", "trial", dv, "group", "nrtrials", "nrchannels"]);
            dat_tbl = stack(dat_tbl, dv, "IndexVariableName", 'dependent', 'NewDataVariableName', 'y');
            dat_tbl= convertvars(dat_tbl,'dependent', 'char');
            dat_tpl = dj.struct.join(table2struct(dat_tbl),key);

            insert(ns.TepochChannel,dat_tpl)

        end
    end

end

function [ecTbl, groups] = find_groups_to_average_(self, key, varargin)

% every char element in varargin is either a custom func_str or one of the
% pre-set averaging options. Every cell array within varargin is associated
% with the preceding char element.

% Get relevant tables
% eTbl = ns.Epoch & key;

n_arg = numel(varargin);
ii = 1;
while ii <= n_arg

    argN = varargin{ii};
    % whether assigned input args
    inp_argsN = {};
    if n_arg > ii && iscell(varargin{ii+1})

        inp_argsN = varargin{ii+1};
        
    end

    % find groups
    switch argN

        case {'trial','channel','condition'}
            
            disp('TBD.');

        otherwise % custom function

            fun = str2func(argN);
            [ecTbl, groups] = fun(self, key, inp_argsN{:});

    end

    % next item
    ii = ii + 1 + ~isempty(inp_argsN);

end


end