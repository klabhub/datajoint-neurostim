function parms = mlspikeDefaults(publication,indicator,fineTuned)
% function parms = mlSpikeDefaults(publication,indicator,fineTuned)
%  Convenience function to setup parms struct with defaults
%  from the literature. 
%
% SEE ALSO sbx.mlspike, Mlcalibration

arguments
    publication (1,1) string
    indicator (1,1) string
    fineTuned (1,1) logical =false
end

parms = tps_mlspikes('par');
switch (publication)
    case "KPOP25"
        % Values that gave good qualtiy for KPOP data
        % (I ran autocalibrate with RUPPRECHT25 on a large set 
        % and then took the average a/tau/)
        switch upper(indicator)
            case "GCAMP6S"
                parms.a = 0.14;  % Amplitude
                parms.tau  =0.97; % Decay tau in s
                % From rupprecht25:
                parms.hill = 1.84;
                parms.ton = 0.07; % Rise tau (t_on) in s
                parms.pnonlin = 0.1;
                parms.drift.parameter = 0.1;
            otherwise 
                error("No %s defaults for %s",indicator,publication)
        end
    case "RUPPRECHT25"
        % Return parameters for Gcamp6s based on Rupprecht et al. 2025.
        switch upper(indicator)
            case "GCAMP6S"
                parms.hill = 1.84;
                parms.a = 0.113;  % Amplitude
                parms.tau  =1.87; % Decay tau in s
                parms.ton = 0.07; % Rise tau (t_on) in s
                parms.pnonlin = 0.1;
                parms.drift.parameter = 0.1;
            case "GCAMP8S"
                parms.hill = 2.2;
                parms.a = 0.576;  % Amplitude
                parms.tau  =0.267; % Decay tau in s
                parms.ton = 0.00472; % Rise tau (t_on) in s
                parms.pnonlin = 0.1;
                parms.drift.parameter = 0.1;
                if fineTuned
                    % Use tau and amplitude as estimated by Rupprecht
                    % with a grid search on ground truth data
                    parms.a    = 1;  % Amplitude
                    parms.tau  = 0.7; % Decay tau in s
                end
            otherwise
                error("No %s defaults for %s",indicator,publication)
        end
    case "DENEUX16"
        switch upper(indicator)
            case "GCAMP6S"
                parms.hill = 1; % Using polynomial pnonlin instead
                parms.a = 0.07;  % Amplitude
                parms.tau  =1.3; % Decay tau in s
                parms.ton = 0.02; % Rise tau (t_on) in s
                parms.pnonlin =  [0.73, -0.05];  % Deneux values for Gcamp6s;
                parms.drift.parameter = 0.1;
            otherwise
                error("No %s defaults for %s",indicator,publication)
        end
    otherwise
        error("Unknown publication %s",publication)
end
end