%{
# Preprocessed epoch data, with epochs stored per trial and channel in the part table ns.EpochChannel
-> ns.C             # The continuous data that were epoched
-> ns.EpochParm     # Parameters used to epoch
-> ns.Dimension     # Dimension that determines the conditions
---
time : blob             # Time in milliseconds relative to the align event (which is defined in EpochParm) [start stop nrSamples]
prep : blob             # Struct with information on preprocessing (.prep) done during epoching
art   : blob            # Struct with information on artifact removal (.art) done during epoching
%}
classdef Epoch < dj.Computed & dj.DJInstance
    properties (GetAccess =public,SetAccess = protected)
        cache table
        cacheQry (1,1) string =""
    end
    properties (Dependent)
        time
        samplingRate
        keySource
    end

    methods %Set/Get
        function v = get.keySource(~)
            % Restricted to Dimensions listed in EpocjParm
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
            v = (proj(ns.C) * proj(ns.EpochParm) * proj(ns.Dimension)) & combinedWhere;
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


    methods (Access = protected)
        function makeTuples(tbl, key)
            %% Determine events to align to
            parmTpl = fetch(ns.EpochParm &key,'prepparms','artparms','align','window','channels');
            conditionTpl = fetch(ns.DimensionCondition&key,'trials');
            trialsInDimension = cat(1,conditionTpl.trials);
            alignTpl = get(ns.Experiment &key,parmTpl.align.plugin,prm=parmTpl.align.event,what=["trialtime" "data" "trial"],trial=trialsInDimension);

            noSuchEvent = isinf(alignTpl.trialtime);
            if any(noSuchEvent)
                fprintf('Removing %d trials in which the %s.%s event did not occur.\n',sum(noSuchEvent),parmTpl.align.plugin,parmTpl.align.event);
                alignTpl.data(noSuchEvent) = [];
                alignTpl.trial(noSuchEvent) =[];
                alignTpl.trialtime(noSuchEvent) =[];
            end
            [trials,ia] = unique(alignTpl.trial,'stable','last');
            if numel(trials) < numel(alignTpl.trial)
                fprintf('The %s event in %s occurs more than once (%d times). Using the last occurrence.\n', parmTpl.align.event,parmTpl.align.plugin,numel(alignTpl.trial) - numel(trials));
            end
            if numel(trialsInDimension) > numel(trials)
                fprintf('%d trials do not have the %s event in %s. Those will not be inclided.\n', numel(trialsInDimension) - numel(trials),parmTpl.align.event,parmTpl.align.plugin);
            end
            startTime = alignTpl.trialtime(ia);
            nrTrials = numel(trials);
            nrConditions = numel(conditionTpl);
            condition = repmat("",nrTrials,1);
            for c=1:nrConditions
                condition(ismember(trials,conditionTpl(c).trials)) = conditionTpl(c).name;
            end

            C = ns.C & key;
            if isempty(parmTpl.channels)
                % Use all channels by default
                parmTpl.channels = C.channels';
            end

            %% Extract aligned segments from ns.C
            tic;
            fprintf("Collecting segmented data from %d channels in ns.CChannel...\n",numel(parmTpl.channels));
            [T,~,channelsWithData] = align(ns.C & key,align=startTime,start=parmTpl.window(1),stop=parmTpl.window(2),trial=trials,channel=parmTpl.channels);
            parmTpl.channels =channelsWithData(:)';
            fprintf("\t Exporting is complete after %s\n",toc);

            %% --- Preprocess epochs ---
            tic;
            fprintf("Preprocessing segmented data...\n");
            [signal,t] = timetableToDouble(T); % [timepoints trials channels ]
            [signal,t,prepResults] = prep.preprocess(signal,seconds(t),parmTpl.prepparms,key);
            [nrSamples,nrTrials,nrChannels] = size(signal); %#ok<ASGLU>
            fprintf("\t Filtering is complete after %s\n",toc);
            %% --- Artifact/Outlier Rejection ---
            tic;
            if isstruct(parmTpl.artparms) || (islogical(parmTpl.artparms) && parmTpl.artparms)
                fprintf("Artifact detection ...\n");
                parmTpl.artparms.epoch_no = trials;
                pv =namedargs2cell(parmTpl.artparms);
                [artResults] = prep.artifactDetection(permute(signal,[2 3 1]),C.samplingRate,pv{:});
                % Remove epochs that were identified as having artifacts
                out = ismember(trials,artResults.all);
                signal(:,out,:) = [];
                trials(out) = [];
                condition(out) = [];
                startTime(out) = [];
                nrTrials =numel(trials);
                fprintf("\t Artifact detection complete after %s\n",toc);
            else
                artResults=struct('dummy',true);
            end

            %% --- Submit to the server ---
            tic;
            fprintf("Submitting to the server\n");
            epoch_tpl = mergestruct(key, ...
                struct(time = [t(1) t(end) numel(t)],...
                prep = prepResults,...
                art =artResults));
            % Insert to Epoch table
            epoch_tpl = makeMymSafe(epoch_tpl);
            insert(tbl, epoch_tpl);

            % Create EpochChannel tuple that contains the data
            signal = reshape(squeeze(num2cell(signal,1)),nrTrials*nrChannels,1);
            trial = num2cell(repmat(trials,nrChannels,1));
            condition = cellstr(repmat(condition,nrChannels,1));
            onset = num2cell(repmat(startTime,nrChannels,1));
            channel  = num2cell(reshape(repmat(parmTpl.channels,nrTrials,1),nrTrials*nrChannels,1));

            if isempty(trial)
                fprintf('No epochs remaining after artifact detection');
            else
                tpl = mergestruct(key, ...
                    struct(signal =signal,...
                    trial = trial, ...
                    onset = onset,...
                    channel = channel,...
                    condition = condition));

                chunkedInsert(ns.EpochChannel, tpl);
                fprintf("\t Submission is complete after %s.\n",toc);
            end
        end
    end
end

