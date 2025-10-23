%{
# Preprocessed epoch data, the actual epochs are stored in ns.EpochChannel
-> ns.C             # The continuous data that were epoched
-> ns.EpochParm     # Parameters used to epoch
-> ns.Dimension     # Dimension that determines the conditions
---
time : blob             # Time in milliseconds relative to the align event (which is defined in EpochParm) [start stop nrSamples]
info : blob             # Info struct with information on preprocessing (.prep) and artifact removal (.art) done during epoching
%}
classdef Epoch < dj.Computed & dj.DJInstance
    properties (Dependent)
        time
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
                parmClauses{p} = sprintf('(ctag="%s" AND dimension = "%s")', ...
                    allParms(p).ctag,allParms(p).dimension);
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


    end

    methods (Access= public)

        function plot (tbl,pv)
            % Plot evoked responses. 
            arguments
                tbl (1,1) ns.Epoch
                pv.averageChannels (1,1) logical =true; % Average over channels
                pv.averageTrials (1,1) logical =true;   % Average over trials
                pv.delta (1,1) string = ""              % Also show the difference using this named condition as the reference
                pv.channel (:,1) double = []            % Select a subset of channels 
            end
            if isempty(pv.channel)
                channelRestrict = "true";
            else
                channelRestrict = struct('channel',num2cell(pv.channel));
            end
            averageDim = [];
            if pv.averageTrials 
                averageDim = [averageDim 2];
            end
            if pv.averageChannels
                averageDim = [averageDim 3];
            end


            for tpl = fetch(tbl,'time')'
                figByName(sprintf("Epoch %s for %s on %s @ %s",tpl.etag,tpl.subject,tpl.session_date,tpl.starttime));
                clf
                hold on
                parms = fetch(ns.EpochParm & tpl,'*');
                t = linspace(tpl.time(1),tpl.time(2),tpl.time(3));
                nrSamples = tpl.time(3);

                conditions = fetch(ns.DimensionCondition & tpl,'name');
                nrConditions = numel(conditions);
                h = gobjects(nrConditions+(pv.delta~="")*(nrConditions-1),1);
                colors ='rgbcmyk';
                for c= 1:nrConditions

                    channelTpl =fetch(ns.EpochChannel &  channelRestrict & tpl & struct('condition',conditions(c).name),'signal','ORDER BY channel');
                    signal = cat(2,channelTpl.signal);
                    channels= [channelTpl.channel];
                    trials = [channelTpl.trial];
                    [~,uChannels] = findgroups(channels);
                    nrChannels = numel(uChannels);
                    [~,uTrials] = findgroups(trials);
                    nrTrials = numel(uTrials);
                    signal = reshape(signal,nrSamples,nrTrials,nrChannels);
                    if isempty(averageDim)
                        se= zeros(size(signal));
                    else
                        se  = std(signal,0,averageDim,"omitmissing")./sqrt(prod(size(signal,averageDim)));
                        signal = mean(signal,averageDim,"omitmissing");
                    end
                    % Reshape to 2D; columsn can be trials (averaged over
                    % channels), channels (averaged over trials) or all
                    % channels and trials (no averaging).
                    signal =reshape(signal,nrSamples,[]);
                    se  =reshape(se,nrSamples,[]);
                    thisH = plot(t,signal,'Color',colors(mod(c-1,numel(colors))+1));
                    h(c) = thisH(1);
                    plot(t,signal+se,':','Color',h(c).Color)
                    plot(t,signal-se,':','Color',h(c).Color)
                    title (tpl.dimension)
                    ylabel 'EP (\muV)'
                    xlabel (sprintf('Time after %s.%s (ms)',parms.align.plugin,parms.align.event));
                end

                if pv.delta ~=""
                    reference = find(pv.delta =={conditions.name});
                    refY = h(reference).YData;
                    others = setdiff(1:nrConditions,reference);
                    deltaLegStr = cell(1,numel(others));
                    for c = 1:numel(others)
                        y = h(others(c)).YData;
                        h(nrConditions+c) = plot(t,y-refY,'Color',colors(mod(nrConditions+c-1,numel(colors))+1));
                        deltaLegStr{c} = [conditions(others(c)).name ' - ' conditions(reference).name ];
                    end
                else
                    deltaLegStr = {};
                end

                cStr = strcat({[tpl.dimension ':']},{conditions.name});
                legend(h,cat(2,cStr,deltaLegStr));
                xlim([min(t) max(t)])
            end

        end
    end
    methods (Access = protected)
        function makeTuples(tbl, key)

            %% Determine events to align to
            parmTpl = fetch(ns.EpochParm &key,'prep','art','align','window','channels');
            conditionTpl = fetch(ns.DimensionCondition&key,'trials');
            trialsInDimension = cat(1,conditionTpl.trials);

            alignTpl = get(ns.Experiment &key,parmTpl.align.plugin,prm=parmTpl.align.event,what=["trialtime" "data" "trial"],trial=trialsInDimension);

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

            if isempty(parmTpl.channels)
                %
                C = ns.C & key;
                parmTpl.channels = C.channels';
            end

            %% Extract aligned segments from ns.C
            tic;
            fprintf("Collecting segmented data from %d channels in ns.CChannel...\n",numel(parmTpl.channels));
            T = align(ns.C & key,align=startTime,start=parmTpl.window(1),stop=parmTpl.window(2),trial=trials,channel=parmTpl.channels);
            fprintf("\t Exporting is complete after %s\n",toc);


            %% --- Preprocess epochs ---
            tic;
            fprintf("Preprocessing segmented data...\n");
            [signal,t] = timetableToDouble(T); % [timepoints trials channels ]
            [signal,t,prepInfo] = ns.prep.preprocess(signal,seconds(t),parmTpl.prep,key);
            [nrSamples,nrTrials,nrChannels] = size(signal); %#ok<ASGLU>
            fprintf("\t Filtering is complete after %s\n",toc);
            %% --- Artifact/Outlier Rejection ---
            tic;
            fprintf("Artifact detection ...\n");
            % samplingRate = 1./mode(diff(time));
            %           [signal,time,artInfo] = ns.prep.artifactDetection(signal,samplingRate,parmTpl.art);
            artInfo =struct('dummy',true);
            fprintf("\t Artifact detection complete after %s\n",toc);

            %% --- Submit to the server ---
            tic;
            fprintf("Submitting to the server\n");
            epoch_tpl = mergestruct(key, ...
                struct(time = [t(1) t(end) numel(t)],...
                info = struct('prep',prepInfo,'art',artInfo)));
            % Insert to Epoch table
            chunkedInsert(tbl, epoch_tpl);

            % Create EpochChannel tuple that contains the actual epoch data
            signal = reshape(squeeze(num2cell(signal,1)),nrTrials*nrChannels,1);
            trial = num2cell(repmat(trials,nrChannels,1));
            condition = cellstr(repmat(condition,nrChannels,1));
            onset = num2cell(repmat(startTime,nrChannels,1));
            channel  = num2cell(reshape(repmat(parmTpl.channels,nrTrials,1),nrTrials*nrChannels,1));

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

