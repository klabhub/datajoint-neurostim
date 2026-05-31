%{
# Labels for independent components 
-> ns.Ica
-> ns.LabelParm
---
q =NULL : longblob  # method-specific data used for labeling (e.g. ICLabel most likely class or ETA-based z-scores)
extra =NULL : longblob # method-specific extra data (e.g. ICLabel full class probabilities or ETA-based event-triggered averages)
%}
%
% After running ICA (see ns.Ica and ns.IcaParm), we want to label
% components for their subsequent use in epoching (ns.Epoch) or to remove
% components during preprocessing.
%
% This table can be filled by different methods labeling the components.
% For instance the eeglab iclabel function will fill q with the highest
% probablity label for each component, and extra with the full probability
% matrix. See ns.Ica/plot() and ns.LabelParm how to set this up for
% populate.
%
% The eta method identifies labels by determining an event triggered average (eta)
% for a (set of) specified events in a plugin and assessing whether that
% event led to a signficant change in the componennt.  This is, for
% instance, a good way to identify blink related components.
% This methods stores a z-score per component; higher z-scores suggest that
% the event really drove the component. (z-score is the difference between
% the component activity at the time of the event compared to a baseline
% period just before).The extra field stores the baseline corrected event
% triggered averages for each component.
% 
% See also ns.LabelParm
%
% BK - May 2026.

classdef Label <  dj.Computed & dj.DJInstance

    properties (Dependent)
        keySource
    end

    methods
        function v = get.keySource(~)
            % Produce (Ica, LabelParm) pairs respecting paradigm restrictions.
            %
            % For each ltag:
            %   - No rows in LabelParmParadigm → applies to all paradigms
            %   - Has rows in LabelParmParadigm → applies only to those paradigms
            %
            % Base query: cross-product of Ica x LabelParm, with Experiment
            % joined to expose the 'paradigm' attribute for filtering.
            allParms = fetch(ns.LabelParm, 'ltag');

            if isempty(allParms)
                v = ns.Ica * proj(ns.Experiment, 'paradigm') * ns.LabelParm & 'FALSE';
                return;
            end

            parmClauses = cell(numel(allParms), 1);
            for p = 1:numel(allParms)
                ltag = allParms(p).ltag;
                paradigms = fetchn(ns.LabelParmParadigm & struct('ltag', ltag), 'name');
                if isempty(paradigms)
                    % Unrestricted: applies to every paradigm
                    parmClauses{p} = sprintf('(ltag="%s")', ltag);
                else
                    % Restricted: only for the listed paradigms
                    paradigmClauses = cellfun(@(pg) sprintf('paradigm="%s"', pg), ...
                        paradigms, 'UniformOutput', false);
                    parmClauses{p} = sprintf('(ltag="%s" AND (%s))', ltag, ...
                        strjoin(paradigmClauses, ' OR '));
                end
            end

            combinedWhere = ['(' strjoin(parmClauses, ' OR ') ')'];
            v = (ns.Ica * proj(ns.Experiment, 'paradigm') * ns.LabelParm) & combinedWhere;
        end
    end


    methods (Access=public)
        function T  = find(tbl,q)
            % Find components matching a query string (for iclabel) or numeric threshold (for eta method).
            arguments
                tbl (1,1) ns.Label
                q (1,:) 
            end
            T = fetchtable(tbl*ns.LabelParm,'*');
            T =addvars(T,cell(height(T),1),'NewVariableNames','components');
            for tpl=1:height(T) 
                if isstruct(T.parms)
                    method = T.parms(tpl).method;
                else
                    method = T.parms{tpl}.method;
                end
                switch upper(method)
                    case 'ICLABEL'
                        if ~iscellstr(q) && ~isstring(q) && ~ischar(q)
                            continue;
                        end
                        comp = find(contains(T.q{tpl}, q,'IgnoreCase',true));
                        if isempty(comp)
                            fprintf('No components matching "%s".\n', q);
                        else
                            fprintf('Components matching "%s": %s\n', q, strjoin(string(comp), ', '));
                        end
                    case 'ETA'
                        if ~isnumeric(q)
                            continue;
                        end
                        comp = find(T.q{tpl} > q);  
                        if isempty(comp)
                            fprintf('No components with z>=%.2f .\n', q);
                        end
                    otherwise 
                        error('No find implemented for labeling method: %s', T.parms{tpl}.method);

                end
                T{tpl,'components'} = {comp};
            end
            

        end
      
    end




    methods (Access=protected)
        function makeTuples(tbl,key)
            % Populate Label for the given key (Ica, LabelParm) tuple.
            % The makeTuples function will be called once for each key in
            % keySource that does not already have a corresponding entry in
            % Label. The key will include the identifying attributes of both
            % ns.Ica and ns.LabelParm, as well as the 'paradigm' attribute
            % from ns.Experiment.
            %
            % This function should compute the labels for the ICA components
            % based on the parameters specified in LabelParm and insert a new
            % tuple into Label with the results.

            parms = fetch1(ns.LabelParm & key,'parms');

            switch upper(parms.method)
                case 'ICLABEL'
                    % Check for eeglab and iclabel function
                    assert(exist('iclabel', 'file') == 2, 'EEGLAB with iclabel function is required for method ''iclabel''.');
                    EEG =  ephys.eeglab.dataset(key,data=key.ctag,itag=key.itag);
                    EEG = iclabel(EEG);
                    extra = EEG.etc.ic_classification.ICLabel.classifications; % The full probability matrix 
                    [~,ix] = max(extra,[],2);
                    q = {EEG.etc.ic_classification.ICLabel.classes(ix)}; % The highest probability class label
                case 'SPEARMAN'
                    EEG =  ephys.eeglab.dataset(key, data=key.ctag,itag=key.itag);
                    c = ns.C & rmfield(key, ["ctag" "filename"]) & struct('ctag', parms.ctag);
                    tpl = fetch(c * ns.CChannel & struct('name', parms.channel), '*');

                    eegSampleTime = polyval(EEG.etc.neurostim.clockParms, EEG.times/1000)';
                    cSampleTime   = linspace(tpl.time(1), tpl.time(2), tpl.time(3))';
                    cSignal = interp1(cSampleTime, tpl.signal, eegSampleTime);

                    % Remove samples where the signal is missing or an outlier
                    cSignal(isnan(cSignal)) = tpl.min;

                    % Spearman correlation: robust to non-Gaussian pupil distribution
                    [q , extra] = corr(EEG.icaact', cSignal, 'Type', 'Spearman');
                        
                case 'ETA'
                    parms.events = string(parms.events);
                    EEG =  ephys.eeglab.dataset(key,data=key.ctag,itag=key.itag);
                    preSamp  = round(parms.window  * EEG.srate / 1000);
                    postSamp = round(parms.window* EEG.srate / 1000);
                    nSamps   = preSamp + postSamp + 1;
                    timesMs  = (-preSamp:postSamp) / EEG.srate * 1000;   % time axis (ms)

                    EEG = ephys.eeglab.addEvents(EEG,parms.plugin,parms.events);
                    if ~any(ismember(parms.events, {EEG.event.type}))
                        warning("No %s events found in this dataset.", strjoin(parms.events));
                        q = [];
                        extra = [];
                        
                    else
                   
                    nComps     = size(EEG.icaact, 1);
                    allEventTypes = erase(lower(string({EEG.event.type})), "'");                        
                    keep = ismember(allEventTypes,parms.events);
                    latencies = round([EEG.event(keep).latency]);                    
                    ok = latencies > preSamp & latencies <= (EEG.pnts - postSamp);
                    latencies = latencies(ok);
                    nEvents = numel(latencies);
                    idx =  ((-preSamp:postSamp) + latencies(:))';                    
                    epochs = reshape(EEG.icaact(:,idx),[nComps nSamps nEvents]);
                     
                    % Baseline window: first few ms of the window
                    baseWin = timesMs >= -parms.window & timesMs < -parms.window+parms.baseline;
                    mu_base = mean(epochs(:, baseWin,:), [2 3]);   
                    sd_base = std(epochs(:, baseWin,:), 0, [2 3]); 
                    mu_resp = mean(epochs(:, ~baseWin,:), [2 3]);  

                    % z-score per component per event type
                    q = (mu_resp - mu_base) ./ sd_base;         

                    % Baseline-corrected ETA: mean over events at each time point [nComps x nSamps]
                    extra = mean(epochs, 3, "omitmissing") - mu_base;
                    end
                otherwise
                    error('Unknown labeling method: %s', parms.method);
            end

            tpl = mergestruct(key,... 
                            struct('q', q, 'extra', extra));
            insert(tbl, tpl)

        end
    end


end



