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

    methods (Access=public)
        function plot(tbl,varargin)
    % Wrapper to call plot on the ns.TEpochChannel table, which is a cache table that contains the actual data.
    % The TEpoch table contains the metadata, but the actual data is in  TEpochChannel.
    channelTbl = ns.TepochChannel & tbl;
    plot(channelTbl,varargin{:});
    
    end
    end
    methods (Access = protected)
        function makeTuples(self, key)
            % Apply a computation/transform to a collection of Epochs and
            % store as Tepoch.
            parms = fetch1(ns.TepochParm & key,'parms');
            fun   = fetch1(ns.TepochParm & key,'fun');
            parms = namedargs2cell(parms);
            % T containts the independent and dependent variables, D maps each independent-variable column to its dependent columns.
            [T,D] = compute( ns.EpochChannel&key,fun,parms{:});
            insertTuples(self,T,D,key)
        end


        function insertTuples(self,T,D,key)
            % Insert the Tepoch and TepochChannel tuples represented by T and D.

            % Extract each dependent variable and insert its Tepoch tuple.
            mapKeys = string(keys(D));
            allDv = string.empty(1,0);
            for iMap = 1:numel(mapKeys)
                idv = mapKeys(iMap);
                dv = string(D(char(idv)));
                allDv = [allDv dv(:)']; %#ok<AGROW>
                for iDv = 1:numel(dv)
                    x = table2cell(T(1,idv)); % Each row should have the same idv
                    x = cat(2,x{:});
                    tpl = dj.struct.join(struct(independent = idv, ...
                        x = x, dependent = char(dv(iDv))),key);
                    insert(self,makeMymSafe(tpl));
                end
            end

            % Collect the information per channel/trial.
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
                T.trial = zeros(height(T),1); % grouped/averaged
            end

            % Build one TepochChannel row per dependent variable and group.
            % Dependent arrays may have different widths, so stack() cannot combine
            % them into one homogeneous table variable.
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




