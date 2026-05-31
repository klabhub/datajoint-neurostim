function   varExplainedPct = icaVarianceExplained(EEG)
            % Helper function to compute variance explained per IC.

            % Channels used for ICA
            X = EEG.data(EEG.icachansind, :);          % [nChan x nTime]

            % ICA activations (sources)
            % Prefer EEG.icaact if already present; otherwise compute from weights/sphere
            if isfield(EEG, 'icaact') && ~isempty(EEG.icaact)
                S = EEG.icaact;                         % [nComp x nTime]
            else
                S = (EEG.icaweights * EEG.icasphere) * X;
            end

            % Mixing matrix
            A = EEG.icawinv;                            % [nChan x nComp]

            % Remove source means for variance calculations
            S0 = S - mean(S, 2);                        % [nComp x nTime]

            % Sensor-space variance contributed by each component:
            % var_i = sum over channels of var(a_i * s_i)
            %       = ||a_i||^2 * var(s_i)
            nTime = size(S0, 2);
            varS = sum(S0.^2, 2) / (nTime - 1);         % [nComp x 1]
            normA2 = sum(A.^2, 1)';                     % [nComp x 1]
            varPerComp = normA2 .* varS;                % [nComp x 1]

            % Total sensor variance of observed data (same channel set)
            X0 = X - mean(X, 2);
            totalVar = sum(sum(X0.^2, 2) / (nTime - 1));

            % Percent variance explained per component
            varExplainedPct = 100 * varPerComp / totalVar;   % [nComp x 1]

        end