%{
# Preprocessed epoch data, with epochs stored per trial and channel in the part table ns.EpochChannel
-> ns.C             # The continuous data that were epoched
-> ns.EpochParm     # Parameters used to epoch
-> ns.Dimension     # Dimension that determines the conditions and selects the trials
---
time : blob             # Time in milliseconds relative to the align event (which is defined in EpochParm) [start stop nrSamples]
prep : blob             # Struct with information on preprocessing (.prepparms) done during epoching
art   : blob            # Struct with information on artifact removal (.artparms) done during epoching
plg   : blob            # Struct with information on epoch removal (.plgparms) based on behavior/plugins done during epoching
%}
classdef Epoch < dj.Computed & dj.DJInstance
   
    properties (Dependent)
        time
        samplingRate
        keySource
    end

    methods %Set/Get
        function v = get.keySource(~)
            % Restricted to Dimensions listed in EpochParm
            % Fetch all TuningParm rows at once
            allParms = fetch(ns.EpochParm, 'etag', 'dimension','ctag');

            if isempty(allParms)
                % No EpochParm rows - return empty result
                v = proj(ns.C) * ns.EpochParm * proj(ns.Dimension) & 'FALSE';
                return;
            end

            % Build a combined WHERE clause for all EpochParm rows using OR
            parmClauses = cell(numel(allParms), 1);

            for p = 1:numel(allParms)
                % Each clause needs both epochtag and dimension
                parmClauses{p} = sprintf('(ctag="%s" AND dimension = "%s" AND etag ="%s")', ...
                    allParms(p).ctag,allParms(p).dimension,allParms(p).etag);
            end

            % Combine all clauses with OR
            combinedWhere = ['(' strjoin(parmClauses, ' OR ') ')'];

            % Apply combined restriction
            % Selecting only those dimensions that have actual conditions.
            v = (proj(ns.C) * proj(ns.EpochParm) * proj(ns.Dimension & ns.DimensionCondition)) & combinedWhere;
        end
        function t =get.time(tbl)
            t = fetchn(tbl, 'time');
            t = cellfun(@(x) linspace(x(1),x(2),x(3))',t,'UniformOutput',false);
            if count(tbl)==1
                t= t{1};
            end
        end
        function v = get.samplingRate(tbl)
            t = fetchn(tbl, 'time');
            v= cellfun(@(x) x(3)./(x(2)-x(1)),t,'UniformOutput',true);
        end        
    end


    methods
        function T = attrition(tbl)
            % Return a table that shows the attrition (i.e. how many trials
            % and channels were discarded per epoch and why)
            T=fetchtable(tbl,'art','plg');  % Get the structs that store the badBy info
            nrEpochs = height(T);
            % Count trials in the dimension for each epoch row
            dimCounts = fetchtable(aggr(tbl, ns.DimensionTrial, 'count(*)->nrInDimension'), 'nrInDimension');
            T = innerjoin(T,dimCounts);

            % Count unique trials in EpochChannel for each epoch row
            epochCounts = fetchtable(aggr(tbl, ns.EpochChannel, 'count(distinct trial)->nrInEpoch'), 'nrInEpoch');
            T = innerjoin(T,epochCounts);


            T = addvars(T,zeros(nrEpochs,1),zeros(nrEpochs,1),'NewVariableNames',{'badByArtifact','badByAlignTime'});

            % For each epoch check why the trials were removed by
            % inspecting the badBy fields of the art and plg structs.
            for cntr = 1:nrEpochs
                art = T.art(cntr);
                if ~isempty(art)
                    T.badByArtifact(cntr)= numel(art.all);
                end
                plg = T.plg(cntr);                
                if ~isempty(plg)
                    for c = string(plg.categories)'
                        c=deblank(c);
                        if ~ismember(c,T.Properties.VariableNames)
                            T= addvars(T,zeros(nrEpochs,1),'NewVariableNames',c);
                        end
                        T.(c)(cntr)= numel(plg.(c));
                    end
                end                
            end

            nrTrialsRemoved = T.nrInDimension-T.nrInEpoch;
            nrBadBy = sum(T{:,startsWith(T.Properties.VariableNames,'badBy')},2);
            isMismatch = nrTrialsRemoved ~=nrBadBy;
            if any(isMismatch)
                fprintf(2,"Trial attrition counts do not match the badBy counts in: \n")
                T(isMismatch,:);
            end
        end
    end

    methods (Access = protected)
        function makeTuples(tbl, key)
            %% Determine events to align and select trials based on plugins
            parmTpl = fetch(ns.EpochParm &key,'prepparms','artparms','plgparms','align','window','channels');
            conditionTpl = fetch(ns.DimensionCondition&key,'trials');
            trialsInDimension = cat(1,conditionTpl.trials);
            alignTpl = get(ns.Experiment &key,parmTpl.align.plugin,prm=parmTpl.align.event,what=["trialtime" "data" "trial"],trial=trialsInDimension);
            if isempty(alignTpl)
                error('This experiment does not use the %s plugin',parmTpl.align.plugin);
            end
            noSuchEvent = isinf(alignTpl.trialtime);
            if any(noSuchEvent)
                fprintf('Removing %d trials in which the %s.%s event did not occur.\n',sum(noSuchEvent),parmTpl.align.plugin,parmTpl.align.event);
                badByAlignTrial = alignTpl.trial(noSuchEvent);
                alignTpl.data(noSuchEvent) = [];
                alignTpl.trial(noSuchEvent) =[];
                alignTpl.trialtime(noSuchEvent) =[];
            else
                badByAlignTrial = [];
            end

            % Select trials based on behavior/plugin parameters
            badByPlg = prep.pluginState(ns.Experiment& key,unique(alignTpl.trial),parmTpl.plgparms);
            outBasedOnPlg  = ismember(alignTpl.trial,badByPlg.all);
            if any(outBasedOnPlg)
                fprintf('Removing %d trials based on plugin parameter selection (%s).\n',sum(outBasedOnPlg ),strjoin(setdiff(fieldnames(parmTpl.plgparms),{'enable'}),'/'));
                alignTpl.data(outBasedOnPlg) = [];
                alignTpl.trial(outBasedOnPlg) =[];
                alignTpl.trialtime(outBasedOnPlg) =[];
            end

            if ~isempty(badByAlignTrial)
                setProperty(badByPlg,'AlignTime',badByAlignTrial);
            end

            %  If an event occurs more than once, use the last.
            %  R2024a and earlier dont allow combining 'stable', 'last'.
            %  This has the same effect
            allTrials = flip(alignTpl.trial);
            allTrialTimes = flip(alignTpl.trialtime);
            [trials,ia] = unique(allTrials,'stable'); 
            if numel(trials) < numel(allTrials)
                fprintf('The %s event in %s occurs more than once (%d times). Using the last occurrence.\n', parmTpl.align.event,parmTpl.align.plugin,numel(alignTpl.trial) - numel(trials));
            end
            startTime = allTrialTimes(ia);
                       
            C = ns.C & key;
            if isempty(parmTpl.channels)
                % Use all channels by default
                parmTpl.channels = C.channels';
            end
            %% Extract aligned segments from ns.C
            tic;
            fprintf("Collecting segmented data from %d channels in ns.CChannel...\n",numel(parmTpl.channels));
            % NOte that this uses the unique trials in which the event
            % occurred and the corresponding align times. These are not
            % stored ascending, but align resorts (and can remove trials
            % too if there are artifacts for instance)
            [T,~,channelsWithData] = align(ns.C & key,align=startTime,start=parmTpl.window(1),stop=parmTpl.window(2),trial=trials,channel=parmTpl.channels);
            
            parmTpl.channels =channelsWithData(:)';
            % Extract the actual trials (in order of signal cols) that have
            % been extracted
            trials = T.Properties.CustomProperties.trials'; 
            startTime = T.Properties.CustomProperties.alignTime';
            
            fprintf("\t Segmenting is complete after %s\n",toc);

            %% --- Preprocess epochs ---
            tic;
            fprintf("Preprocessing segmented data...\n");
            [signal,t] = timetableToDouble(T); % [timepoints trials channels ]
            [signal,t,prepResults] = prep.preprocess(signal,seconds(t),parmTpl.prepparms,key);
            [nrSamples,nrTrials,nrChannels] = size(signal); %#ok<ASGLU>
            fprintf("\t Filtering is complete after %s\n",toc);
            %% --- Artifact/Outlier Rejection ---
            tic;
            fprintf("Artifact detection ...\n");
            pv =namedargs2cell(parmTpl.artparms);
            [badByArt] = prep.artifactDetection(permute(signal,[2 3 1]),C.samplingRate,'epoch_no',trials,pv{:});
            % Remove epochs that were identified as having artifacts
            out = ismember(trials,badByArt.all);
            signal(:,out,:) = [];
            trials(out) = [];
            startTime(out) = [];
            nrTrials =numel(trials);

            fprintf("\t Artifact detection complete after %s\n",toc);

            %% --- Submit to the server ---
            tic;
            fprintf("Submitting epochs to the server\n");
            epoch_tpl = mergestruct(key, ...
                struct(time = [t(1) t(end) numel(t)],...
                prep = prepResults,...
                art =badByArt, ...
                plg = badByPlg));
            % Insert to Epoch table
            epoch_tpl = makeMymSafe(epoch_tpl);
            insert(tbl, epoch_tpl);

            % Create EpochChannel tuple that contains the data
            signal = reshape(squeeze(num2cell(signal,1)),nrTrials*nrChannels,1);
            trial = num2cell(repmat(trials,nrChannels,1));
            onset = num2cell(repmat(startTime,nrChannels,1));
            channel  = num2cell(reshape(repmat(parmTpl.channels,nrTrials,1),nrTrials*nrChannels,1));

            if isempty(trial)
                fprintf('No epochs remaining after artifact detection');
            else
                tpl = mergestruct(key, ...
                    struct(signal =signal,...
                    trial = trial, ...
                    onset = onset,...
                    channel = channel));

                chunkedInsert(ns.EpochChannel, tpl);
                fprintf("\t Submission is complete after %s.\n",toc);
            end
        end
    end
end

