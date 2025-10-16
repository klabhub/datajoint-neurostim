%{
# ns.Evoked is average ERPs per channel and dimension
# by default, grouped by dimension conditions
# if a group_fun provided in ns.EvokedParm, it is used to further split
# into groups.
# Part tables, ns.EvokedChannel contain actual signal
-> ns.C
-> ns.DimensionCondition
-> ns.EvokedParm
group : varchar(32)
---
n_epoch: int
%}

classdef Evoked < dj.Computed & dj.DJInstance

    properties

        keySource
        trials_ dj.DJProperty = dj.DJProperty.empty

    end

    properties (Dependent)

        trials

    end


    methods

        function v = get.keySource(varargin)

            % only those events with some clean epochs
            v = (ns.C * ns.DimensionCondition * ns.EvokedParm) & (ns.Epoch & 'flag=""');

        end


    end


    methods (Access = protected)

        function makeTuples(evTbl, key)

            fprintf("\tPreparing to fetch data.\n")
            % --- Get assoc. tables ---
            cTbl = ns.C & key & ns.Epoch; % '& ns.Epoch' to make sure there are epochs assoc.
            c_n_row = count(cTbl);
            assert(c_n_row == 1, "There must be only one C table entry associated with the key, found %d.", c_n_row);
            
            evp_tpl = fetch(ns.EvokedParm & key, '*');

            % ch_qry for later
            ch_qry = sprintf('channel not in (%s)', join(string(cTbl.artifacts.all),','));

            % trial numbers are stored in dimTbl
            dimTbl = ns.DimensionCondition & key;
            dim_tpl = fetch(dimTbl, '*');
            trl_qry = sprintf('trial in (%s)',join(string(dim_tpl.trials),','));

            % epoch table
            eTbl = ns.Epoch & key & trl_qry & evp_tpl.epoch_query; %by default 'flag =""': only clean epochs
            if count(eTbl)==0

                fprintf("\t No clean epochs were found. Skipping...");
                return;

            end
            

            if isempty(evp_tpl.group_fun)

                groups = {'', eTbl};
            else

                % group_fun must accept inputs (eTbl, key) and return a
                % struct scalar whose fields are the split eTbl instances, and whose
                % fieldnames are the group label. Transform this to a cell
                % vector {fieldname1, fieldvalue1,...}
                groups = feval(evp_tpl.group_fun, eTbl, key);
                if isempty(groups)
                    fprintf("\t No clean epochs were found matching with the group description. Skipping...\n");
                    return;
                end
                groups = namedargs2cell(groups);

            end

            % create a tpl per group
            n_gru = numel(groups)/2;
            evk_tpl = repmat(key, n_gru, 1);
            evkch_tpl = cell(1,n_gru);

            for iGru = 1:2:2*n_gru

                t_fetch = tic;
                gru_no = ceil(iGru/2);
                fprintf("\tFetching signal for Group %d/%d and preparing submission...\n", gru_no, n_gru);
                if count(groups{iGru+1})==1;
                    fprintf("\t\t No clean epochs were found. Skipping the group...\n");
                    continue;
                end
                evk_tpl(gru_no).group = groups{iGru};
                evk_tpl(gru_no).n_epoch = count(groups{iGru+1});

                % create EvokedChannel tuple
                ecTbl = ns.EpochChannel & groups{iGru+1} & ch_qry;
                if count(ecTbl)
                    ec_tbl = fetchtable(ecTbl, 'signal');
                    toc(t_fetch);
                    % average by channel
                    [gru_idx, ch_id] = findgroups(ec_tbl.channel);
    
                    evkch_tpl{gru_no} = mergestruct(repmat(rmfield(evk_tpl(gru_no), 'n_epoch'), numel(ch_id), 1),...
                        struct(channel = num2cell(ch_id), ...
                        signal =  mat2cell( ...
                        splitapply(@(x) mean(x,1, 'omitmissing'), ec_tbl.signal, gru_idx), ones(numel(ch_id),1) ...
                        )));
                else
                    fprintf("\t\t No clean epochs associated with the group were found.\n\n");
                end

            end

            evkch_tpl = cat(1, evkch_tpl{:});

            % --- Insert Tables ---
            t_sub = tic;
            fprintf('\t Submitting to the server...\n')
            insert(evTbl, evk_tpl);
            chunkedInsert(ns.EvokedChannel, evkch_tpl);
            toc(t_sub);

        end

    end

   

end

