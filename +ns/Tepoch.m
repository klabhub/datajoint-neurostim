%{
# Transformed Epoch - a computation applied to a (set of) epochs.
-> ns.Epoch                     # The epochs that were transformed
-> ns.TepochParm                # Parameters used for the transformation
dependent       : varchar(64)   # The name of the dependent variable(s)
---
x : blob                        # The values of the independent variable
independent     : varchar(64)   #  The name of the independent variable(s), if multiple, concatenated with ':' in order
%}
classdef Tepoch < dj.Computed & dj.DJInstance

    properties (SetAccess = protected)
        keySource 
    end

    methods 
        function v = get.keySource(self)
            v  = ns.TepochParm * (ns.Epoch & ns.EpochChannel);
        end
    end
    methods (Access = protected)
        function makeTuples(self, key)
            % Apply a computation/transform to a collection of Epochs and
            % store as Tepoch.
            parms = fetch1(ns.TepochParm & key,'parms');                 
            fun   = fetch1(ns.TepochParm & key,'fun');
            parms = namedargs2cell(parms);
            [T,dv,idv] = compute( ns.EpochChannel&key,fun,parms{:});

            
            %% Insert in the table           
            x = table2cell(T(1,idv));
            x = cat(2,x{:});
            tpl = dj.struct.join(struct(independent = strjoin(idv,':'), ...
                x = x, dependent = cellstr(dv(:))),...%, groups=avg_groups),
                key);
            
            insert(self,makeMymSafe(tpl));            

            % dat_tbl = T(:,["channel", "trial", dv, "group", "nrtrials", "nrchannels"]);
            varnames = intersect(["channel" "trial" "nrchannels" "nrtrials" dv "group"],T.Properties.VariableNames);
            dat_tbl = T(:,varnames);
            if ismember("channel",dat_tbl.Properties.VariableNames)
                dat_tbl.nrchannels = ones(height(dat_tbl),1);
            else
                dat_tbl.channel = zeros(height(dat_tbl),1); % Must be grouped/averaged
            end
            if ismember("trial",dat_tbl.Properties.VariableNames)
                dat_tbl.nrtrials = ones(height(dat_tbl),1);
            else
                dat_tbl.trial = zeros(height(dat_tbl),1); % Grouped/Averaged
            end

            dat_tbl = stack(dat_tbl, dv, "IndexVariableName", 'dependent', 'NewDataVariableName', 'signal');
            dat_tpl = dj.struct.join(table2struct(dat_tbl),key);            
            chunkedInsert(ns.TepochChannel,makeMymSafe(dat_tpl))

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

        case '' % no averaging

            ecTbl = ns.EpochChannel & (ns.Epoch & key);
            groups = struct(name = '', id = 1);

        otherwise % custom function

            fun = str2func(argN);
            [ecTbl, groups] = fun(self, key, inp_argsN{:});

    end

    % next item
    ii = ii + 1 + ~isempty(inp_argsN);

end


end
