classdef EpochRestricter

    % A restrictor class used to subset RetrievedEpochs instances
    % restrictor = EpochRestricter(...)
    % restricted_epochs = restrictor.restrict(epochs)

    properties

        channels
        time_window
        frequency_window
        trials
        conditions struct

    end

    methods
        
        function obj = EpochRestricter(pv)

            arguments

                pv.channels = {}
                pv.time_window = {}
                pv.frequency_window = {}
                pv.trials = {}
                pv.conditions = struct()
            
            end

            obj.channels = pv.channels;
            obj.time_window = pv.time_window;
            obj.frequency_window = pv.frequency_window;
            obj.trials = pv.trials;
            obj.conditions = pv.conditions;

        end

        function ep = restrict(obj, ep)

            arguments
                obj (1,1) ns.EpochRestricter
                ep (1,1) ns.RetrievedEpochs
            end

            % Determine non-empty restricter variables
            var_names = ep.data.Properties.VariableNames;

            % Condition vars

            if ~isempty(obj.conditions)

                conds = string(fieldnames(obj.conditions));

                isCondAVar = ismember(conds, var_names);
                assert(all(isCondAVar), ...
                    "No variable named '%s' is in epoch data.", ...
                    join(conds(~isCondAVar), "', '"));

                isCond = ones(1,height(ep.data)) == 1;                
                for iCond = 1:length(conds)

                    condN = conds(iCond);
                    valN = obj.conditions.(condN);

                    isCond = isCond && ismember(ep.data.(condN), valN);

                end

                ep.select_rows(isCond);

            end

            % Data indices
            idx_selector_fq = repmat({':'}, ep.n_rows, max(struct2array(ep.dimensions_)));
            idx_selector_t = repmat({':'}, ep.n_rows, max(struct2array(ep.dimensions_)));

            % Time indices

            if ~isempty(obj.time_window)

                assert(~isempty(ep.temporal_data_columns), "No temporal data is present to subset.")
                isTInWin = gen.ifwithin(ep.timepoints, obj.time_window);
                assert(any(isTInWin), "Time window is out of epoch range.");

            else

                isTInWin = ones(1, ep.epoch_length) == 1;

            end
            idx_selector_t(:,end) = repmat({isTInWin}, 1, ep.n_rows);

            % Frequency indices

            if ~isempty(obj.frequency_window)

                assert(~isempty(ep.spectral_data_columns), "No spectral data is present to subset.")
                isFqInWin = gen.ifwithin(ep.frequencies, obj.frequency_window);
                assert(any(isFqInWin), "Frequency window is out of range.");

            else

                isFqInWin = ones(size(ep.frequencies)) == 1;

            end
            idx_selector_fq(:,end) = repmat({isFqInWin}, 1, ep.n_rows);

            

            % Trial indices
            % ADD FEATURE TO EXCEPT FUNCTION HANDLES, CELL ARRAYS WITH
            % FUNCTION HANDLES, (N_ROWS OR 1, 2) CELL ARRAY WITH FUNCTION HANDLES & ARGS
            if ~isempty(obj.trials)

                if isvector(obj.trials)

                    trls = repmat({obj.trials}, 1, ep.n_rows);

                elseif iscell(obj.trials) && (obj.trials) == 1

                    trls = repmat(obj.trials, 1, ep.n_rows);

                else
                    trls = obj.trials;
                end
                assert(~isnan(ep.dimensions_.trials), "Epoch data table does not contain trials.");
                
                idx_selector_t(:, ep.dimensions_.trials) = cellfun(@(x,y) (isempty(x) & y) | ismember(y, x), trls', ep.data.trials, 'UniformOutput', false);
                idx_selector_fq(:, ep.dimensions_.trials) = idx_selector_t(:, ep.dimensions_.trials);

            end

            
            assert(isempty(obj.channels) || ~isnan(ep.dimensions_.channels), "Epoch data table does not contain channels.");
            % Channel indices
            if ~isempty(obj.channels)

                if isvector(obj.channels)

                    ch = repmat({obj.channels}, 1, ep.n_rows);

                elseif iscell(obj.channels) && (obj.channels) == 1

                    ch = repmat(obj.channels, 1, ep.n_rows);

                else
                    ch = obj.channels;
                end
                assert(~isnan(ep.dimensions_.channels), "Epoch data table does not contain channels.");
                idx_selector_t(:, ep.dimensions_.channels) = cellfun(@(x,y) (isempty(x) & y) | ismember(y, x), ch', ep.data.channels, 'UniformOutput', false);
                idx_selector_fq(:, ep.dimensions_.channels) = idx_selector_t(:, ep.dimensions_.channels);

            end           


            % Implement
            for iVar = 1:length(ep.temporal_data_columns)

                varN = ep.temporal_data_columns(iVar);
                datN = cellfun(@(x, idx) x(idx_selector_t{idx,:}), ...
                    ep.data.(varN), num2cell(1:ep.n_rows)', ...
                    'UniformOutput', false);
                ep = ep.replace(varN, datN);

            end

            for iVar = 1:length(ep.spectral_data_columns)

                varN = ep.spectral_data_columns(iVar);
                datN = cellfun(@(x, idx) x(idx_selector_fq{idx,:}), ...
                    ep.data.(varN), num2cell(1:ep.n_rows)', ...
                    'UniformOutput', false);
                ep = ep.replace(varN, datN);

            end

        end

    end

end