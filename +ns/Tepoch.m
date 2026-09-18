%{
# Transformed Epoch - a computation applied to a (set of) epochs.
-> ns.Epoch                     # The epochs that were transformed
-> ns.TepochParm                # Parameters used for the transformation
dependent       : varchar(64)   # The name of the dependent variable(s)
---
x : blob                        # The values of the independent variable
independent     : varchar(64)   #  The name of the independent variable(s), if multiple, concatenated with ':' in order
%}
%
% See ns.cache for a list of computations , or how to add your own
% computation (defined in a ns.TepochParm)
classdef Tepoch < dj.Computed & dj.DJInstance

    properties (SetAccess = protected)
        keySource 
    end

    methods 
        function v = get.keySource(~)
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
            [T,D] = compute( ns.EpochChannel&key,fun,parms{:});
            % D maps each independent-variable column to its dependent columns.
            % Extract each dv here and insert.
            mapKeys = string(keys(D));
            allDv = string.empty(1,0);
            for iMap = 1:numel(mapKeys)
                idv = mapKeys(iMap);
                dv = string(D(char(mapKeys(iMap))));
                allDv = [allDv dv(:)']; %#ok<AGROW>
                for iDv = 1:numel(dv)
                    x = table2cell(T(1,idv));
                    x = cat(2,x{:});
                    tpl = dj.struct.join(struct(independent = idv, ...
                        x = x, dependent = char(dv(iDv))),key);
                    insert(self,makeMymSafe(tpl));
                end 
            end
            
            % Collect the information per Channel/Trial
            varnames = intersect(["channel" "trial" "nrchannels" "nrtrials" allDv "group"],T.Properties.VariableNames);
            T = T(:,varnames);
            if ismember("channel",T.Properties.VariableNames)
                T.nrchannels = ones(height(T),1);
            else
                T.channel = zeros(height(T),1); % grouped/averaged
            end
            if ismember("trial",T.Properties.VariableNames)
                T.nrtrials = ones(height(T),1);
            else
                T.trial = zeros(height(T),1); % Grouped/Averaged
            end

            % Build one TepochChannel row per dependent variable and group.
            % Dependent arrays may have different widths, so stack() cannot
            % combine them into one homogeneous table variable.
            channelRows = cell(numel(allDv),1);
            for iDv = 1:numel(allDv)
                dv = allDv(iDv);
                values = T.(dv);
                if ~iscell(values)
                    values = num2cell(values,2);
                end
                channelRow = removevars(T,allDv);
                channelRow.dependent = repmat(dv,height(T),1);
                channelRow.signal = values;
                channelRows{iDv} = channelRow;
            end
            T = vertcat(channelRows{:});
            tpl = dj.struct.join(table2struct(T),key);
            chunkedInsert(ns.TepochChannel,makeMymSafe(tpl))

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
