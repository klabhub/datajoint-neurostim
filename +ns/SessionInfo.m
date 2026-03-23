% -------------------------------------------------------------------------
% ns.SessionInfo
% -------------------------------------------------------------------------
%{
# ns.SessionInfo: Computes session ordering per subject and session
-> ns.Subject
-> ns.SessionInfoParms
---
n_sessions: int
n_days: int
%}
classdef SessionInfo < dj.Computed & dj.DJInstance
    
    properties
        % Restrict population to combinations of Subject and Parms that have an entry in ns.C
        keySource = (ns.Subject * ns.SessionInfoParms) & ns.C;
    end
    
    properties (Dependent)
        % Returns master entries that no longer have associated part entries
        orphaned
    end
    
    methods
        function val = get.orphaned(obj)
            val = obj - ns.SessionInfoData;
        end
    end
    
    methods (Access=protected)
        function makeTuples(obj, key)
            % Fetch parameters to ensure we have the target paradigm
            p = fetch(ns.SessionInfoParms & key, '*');
            
            % Fetch primary keys of the specific sessions matching the paradigm. 
            % Because paradigm is a secondary attribute, extracting the keys 
            % explicitly forces the intersection on shared primary keys.
            valid_sessions = fetch(ns.Experiment & sprintf('paradigm="%s"', p.paradigm));
            
            % Define the subset of ns.C for THIS subject matching the parameters AND valid sessions
            cTbl = ns.C & key & valid_sessions;
            % Fetch necessary primary keys for time sorting and insertion
            c_data = fetch(cTbl, 'session_date', 'starttime', 'filename');
            
            if isempty(c_data)
                return;
            end
            
            % Compute datetime for sorting
            date_strs = {c_data.session_date}';
            time_strs = {c_data.starttime}';
            timestamps = datetime(strcat(date_strs, {' '}, time_strs));
            
            % Sort chronologically
            [sorted_timestamps, sort_idx] = sort(timestamps);
            c_data = c_data(sort_idx);
            
            % Calculate session_day (chronological day index)
            dates_only = dateshift(sorted_timestamps, 'start', 'day');
            [~, ~, day_indices] = unique(dates_only, 'sorted');
            
            s_tpl = key;
            s_tpl.n_sessions = length(c_data);
            s_tpl.n_days = max(day_indices);

            % 1. Insert Master tuple
            insert(obj, s_tpl);
            
            
            % 3. Prepare and insert Data part tuples
            data_tuples = repmat(key, length(c_data), 1);
            
            for i = 1:length(c_data)
                % Map specific session primary keys
                data_tuples(i).session_date = c_data(i).session_date;
                data_tuples(i).starttime = c_data(i).starttime;
                data_tuples(i).filename = c_data(i).filename;
                
                % Compute new attributes
                data_tuples(i).session_no = i;
                data_tuples(i).session_day = day_indices(i);
                
                if i == 1
                    data_tuples(i).time_since_prev = '00:00:00';
                else
                    dur_seconds = seconds(sorted_timestamps(i) - sorted_timestamps(i-1));
                    h = floor(dur_seconds / 3600);
                    m = floor(mod(dur_seconds, 3600) / 60);
                    s = floor(mod(dur_seconds, 60));
                    data_tuples(i).time_since_prev = sprintf('%02d:%02d:%02d', h, m, s);
                end
            end
            
            % Insert into Part table
            insert(ns.SessionInfoData, data_tuples);
        end
    end
end
