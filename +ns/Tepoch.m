%{
# Transformed Epoch - a computation applied to a (set of) epochs.
-> ns.Epoch         # The epochs that were transformed
-> ns.TepochParm    # Parameters used for the transformation
dependent : varchar(64)  # The name of the dependent variable(s) - vector of strings 
---
x : blob            # The values of the independent variable
independent : varchar(64)  #  The name of the independent variable - vector of strings
%}
classdef Tepoch < dj.Computed & dj.DJInstance

    methods (Access = protected)
        function makeTuples(tbl, key)
            % Apply a computation/transform to a collection of Epochs and
            % store as Tepoch.
            parms = fetch(ns.TepochParm & key,'*');

            % Restrict the epoch channels and trials if requested in the
            % parms
            ecTbl = (ns.EpochChannel & key);
            if ~isempty(parms.channels)
                ecTbl = ecTbl & struct('channel',num2cell(parms.channels));
            end
            if ~isempty(parms.trials)
                ecTbl = ecTbl & struct('trial',num2cell(parms.trial));
            end
            if isempty(parms.window)
                % Use the same window as the epoch
                parms.window = fetch(ns.EpochParm & key,'window');
            end
            % Compute - uses the ns.cache/compute function  (abstract
            % superclass).
            [T,dv,idv] = compute(ecTbl,parms.fun,average=parms.average,timeWindow=parms.window);

            % Insert in the table
            %%
            x = table2cell(T(1,idv));
            x = cat(2,x{:});
            tpl = dj.struct.join(struct(independent = strjoin(idv,':'), x = x, dependent = cellstr(dv(:))),key);

            % tpl = key;
            % tpl.x = table2cell(T(1,idv)); % Store the first row only-all others are the same (see fill)
            % tpl.dependent = dv;
            % tpl.independent = strjoin(idv,':');
            insert(tbl,tpl);


            if ismember("trial",parms.average)
                % Trials were averaged out (per-condition average).
                % Use a fake trial number that corresponds to the first
                % trial in each condition
                dimTrials = fetchtable(proj(ns.DimensionTrial & key,'name->condition'));
                dimTrials = innerjoin(T,dimTrials);
                G= groupsummary(dimTrials,"condition",{@min,@numel},"trial");
                G=renamevars(G,["fun1_trial" "fun2_trial"],["trial" "count"]);
                T = innerjoin(T,G,"Keys","condition","RightVariables",["trial" "count"]);
                nrTrials = T.count;
            else
                nrTrials = ones(height(T),1);
            end

            if ismember("channel",parms.average)
                % Channel was averaged out, replace by 0
                T = addvars(T,zeros(height(T),1),'NewVariableNames','channel');
                nrChannels = numel(unique([fetch(ecTbl,'channel').channel]))*ones(height(T),1);

            else
                nrChannels =ones(height(T),1);
            end


            dat_tbl = T(:,["channel", "trial", dv]);
            dat_tbl.nrtrials = nrTrials;
            dat_tbl.nrchannels = nrChannels;
            dat_tbl = stack(dat_tbl, dv, "IndexVariableName", 'dependent', 'NewDataVariableName', 'signal');
            dat_tbl= convertvars(dat_tbl,'dependent', 'char');
            dat_tpl = dj.struct.join(table2struct(dat_tbl),key);

            insert(ns.TepochChannel,dat_tpl)

        end
    end

end