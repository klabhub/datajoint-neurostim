%{

# Epoching parameters to segment trial data. All times are relative to the plugin time

etag : varchar (32) # unique tag
paradigm : varchar(32)
dimension : varchar(32) # Condition from the dimension table
---
epoch_win : BLOB
plugin_name : varchar(32) # The plugin to be timelocked to from DimensionCondition table
pv : BLOB # structure array containing the custom parameters

%}
%
% parms
%   parms.epoch_win double (1,2) # The relative start and end time of the epochs with respect to the plugin time
%   parms.baseline bool (1,1) # whether to baseline or not
%   parms.baseline_win double (1,2) # The baseline window (1, 1)
%
% MOz Feb, 2025
classdef EpochParameter < dj.Lookup

    properties (Dependent)

        etag
        paradigm
        dimension
        plugin_name
        epoch_win
        baseline
        baseline_win
        detrend
        nrows

    end

    methods

        function n = get.nrows(eTbl)

            n = count(eTbl);

        end

        function e = get.etag(eTbl)
            
            tbl = fetch(eTbl, 'etag');

            e = string(vertcat(tbl.etag));
            
        end

        function e = get.paradigm(eTbl)
            
            tbl = fetch(eTbl, 'paradigm');

            e = string(vertcat(tbl.paradigm));
            
        end

        function e = get.dimension(eTbl)
            
            tbl = fetch(eTbl, 'dimension');

            e = string(vertcat(tbl.dimension));
            
        end

        function e = get.plugin_name(eTbl)
            
            tbl = fetch(eTbl, 'plugin_name');

            e = string(vertcat(tbl.plugin_name));
            
        end

        function win = get.epoch_win(eTbl)

            tbl = fetch(eTbl, 'epoch_win');

            win = vertcat(tbl.epoch_win);
            

        end

        function win = get.baseline_win(eTbl)

            parms = fetch(eTbl, 'pv');

            win = eTbl.populate_props_(nan(eTbl.nrows,2), {parms.pv}, 'baseline_win');

        end

        function i = get.baseline(eTbl)
            
            parms = fetch(eTbl, 'pv');
            i = eTbl.populate_props_(zeros(eTbl.nrows,1)==1, {parms.pv}, 'baseline');

        end

        function i = get.detrend(eTbl)
            
            parms = fetch(eTbl, 'pv');

            i = eTbl.populate_props_(zeros(eTbl.nrows,1)==1, {parms.pv}, 'detrend');

        end
        
    end

    methods (Access = private, Static)

        function op = populate_props_(op, entity, prop)

            nrow = size(op,1);
            ncol = size(op,2);

            idx = repmat({':'}, 1, ncol);



            for i = 1:nrow

                if iscell(entity(i))
                    entityN = entity{i};
                elseif isstruct(entity(i))
                    entityN = entity(i);
                end

                if isfield(entityN, prop) || isprop(entityN(i), prop)

                    idx{1} = i;
                    op(idx{:}) = entityN.(prop);

                end                

            end

        end

    end

end