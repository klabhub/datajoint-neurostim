%{
# Stores Independent Component Analysis results for a row in ns.C
-> ns.C
-> ns.IcaParm
---
nrcomponents :  int  # Number of components
winverse : longblob  # Inverse ICA weights (mixing matrix)
sphere   : longblob  # Sphering matrix
weights  : longblob  # ICA weights (unmixing)
channels : longblob  # Channels used for ICA
variance : blob      # Variance explained per component.
%}
%
% Currently implemented with EEGLAB
%
% See also ns.IcaParm, ephys.eeglab.preprocess
%
% BK - May 2026
classdef Ica < dj.Computed
    properties (Dependent)
        keySource
    end

    methods
        function v = get.keySource(~)
            % Restricted to ns.C tuples with the ctag specified in IcaParm
            v = (ns.C * ns.CParm) * proj(ns.IcaParm,'itag','ctag');
        end
    end

    methods (Access=public)
        function plot(tbl,pv)
            arguments
                tbl (1,1) ns.Ica
                pv.comp (1,:) double {mustBeInteger,mustBePositive} = 1:4
                pv.layout (1,2) double {mustBeInteger,mustBePositive}= [4 4]
            end
            tpl = fetch(tbl,'*');

            assert(~isempty(tpl),'No rows in ns.Ica to plot.');
            for i = 1:numel(tpl)
                figure;
                this = tpl(i);
                nComp = size(this.winverse,2);
                if isscalar(pv.comp)
                    comps = 1:min(pv.comp,nComp);
                else
                    comps = pv.comp(pv.comp<=nComp);
                end
                assert(~isempty(comps), 'No requested components exist in this row (nComp=%d).', nComp);

                chanlocs = [fetch(ns.CChannel & this,'channelinfo').channelinfo];

                %   chanlocs = sanitizeChanlocs(chanlocs);

                locFile = [tempname '.loc'];
                writelocs(chanlocs,locFile,'filetype','loc','header','off');
                pageTitle = sprintf('%s @ %s %s | %s/%s', ...
                    this.subject,this.session_date,this.starttime,this.ctag,this.itag);
                compmap(this.winverse,locFile,comps,pageTitle,pv.layout);

                if exist(locFile,'file')
                    delete(locFile);
                end
            end
        end
    end
    methods (Access=protected)
        function makeTuples(tbl,key)

            icaParms = fetch1(ns.IcaParm & key,'parms');
            % Get an EEG struct from ns.C using the ctag specified in IcaParm
            EEG = ephys.eeglab.dataset(key,data=key.ctag);
            % Use the parms from IcaParm to run ICA on the EEG struct, which will add the ICA fields to the struct. The parms should be a struct with a field 'eeglab' that contains a struct of parameters to pass to pop_runica. For example:
            parms.eeglab.ica =  icaParms; % Match the field name expected by preprocess.m
            EEG = ephys.eeglab.preprocess(EEG,parms);

            varExplained = ns.Ica.varianceExplained(EEG);


            tpl = key;
            tpl.nrcomponents = size(EEG.icaweights,1);
            tpl.winverse = EEG.icawinv;
            tpl.sphere = EEG.icasphere;
            tpl.weights = EEG.icaweights;
            tpl.channels = EEG.icachansind;
            tpl.variance = varExplained;
            insert(tbl,tpl);
        end
    end

    methods (Static)
        function   varExplainedPct = varianceExplained(EEG)
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
    end
end
