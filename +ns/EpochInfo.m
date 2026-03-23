% -------------------------------------------------------------------------
% ns.EpochInfo
% -------------------------------------------------------------------------
%{
# ns.EpochInfo: Computes the true chronological order of epochs
-> ns.C
-> ns.EpochInfoParm
%}
classdef EpochInfo < dj.Computed
    properties
        % Restrict population to combinations      
        keySource = ns.EpochInfoParm * proj(ns.C & proj(ns.Epoch & proj(ns.EpochInfoParm,'etag')));
    end
    
    properties (Dependent)
        % Returns a MATLAB table with n_trl, n_valid_trl, and pc_valid_trials
        trial_counts
    end
    
    methods
        function res = get.trial_counts(obj)
            % Aggregate total trials
            q_total = aggr(obj, ns.EpochInfoData, 'count(*)->n_trl');
            
            % Aggregate valid trials
            q_valid = aggr(q_total, ns.EpochInfoData & 'trial_order IS NOT NULL', 'count(*)->n_valid_trl');
            
            % Fetch and convert to MATLAB table
            res = proj(q_valid*q_total,'*', '100*n_valid_trl / n_trl -> pc_valid_trl');
        end
    end
    
    methods (Access=protected)
        function makeTuples(obj, key)
            p = fetch(ns.EpochInfoParm & key, '*');
            cTbl = ns.C & (ns.Epoch & proj(ns.EpochInfoParm & key,'etag'));
            eTbl = ns.Epoch & proj(cTbl) & key & sprintf('etag="%s"', p.etag);
            
            if count(eTbl) == 0
                return;
            end
            
            % Restrict to valid trials
            if isempty(p.exclusion_query) || strlength(strtrim(p.exclusion_query)) == 0
                eTbl_vld = eTbl;
            else
                eTbl_vld = eTbl & p.exclusion_query;
            end
            
            % Fetch all trials for the part table, and valid trials for ordering
            all_ep = fetch(eTbl, 'trial');
            vld_ep = fetch(eTbl_vld, 'trial');
            
            % Sort valid trials by onset to establish sequential order
            [~, sort_idx] = sort([vld_ep.trial]);
            vld_ep = vld_ep(sort_idx);
            valid_trial_nums = [vld_ep.trial];
            
            % Insert Master tuple
            insert(obj, key);
            
            % Prepare and insert Part tuples
            part_tuples = repmat(key, length(all_ep), 1);
            for i = 1:length(all_ep)
                t_num = all_ep(i).trial;
                part_tuples(i).trial = t_num;
                
                % Find sequence order if valid, else NaN (DataJoint converts to NULL)
                idx = find(valid_trial_nums == t_num, 1);
                if ~isempty(idx)
                    part_tuples(i).trial_order = idx;
                else
                    part_tuples(i).trial_order = NaN;
                end
            end
            part_tuples = dj.struct.join(fetch(proj(eTbl)), part_tuples);
            insert(ns.EpochInfoData, part_tuples);
        end
    end
end
