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


    end

    methods (Access= public)

        function plot (tbl,pv)
            % Plot evoked responses for all rows in the table.
            % Set the 'average' input to select which aspects to average
            % over. By default, trials and channels are averaged. 
            % 
            arguments
                tbl (1,1) ns.Epoch {mustHaveRows}
                pv.delta (1,1) string = ""              % Show the difference using this named condition as the reference
                pv.channel (:,1) double = []            % Select a subset of channels
                pv.average (1,:) string {mustBeMember(pv.average,["starttime" "condition" "trial" "channel" "subject" "session_date"])} = ["trial" "channel"]  % Average over these dimensions
                pv.tilesPerPage (1,1) double = 6        % Select how many tiles per page.
            end
            if isempty(pv.channel)
                channelRestrict = "true"; % No restriction
            else
                channelRestrict = struct('channel',num2cell(pv.channel));
            end
            % Fetch the data
            T =fetchtable(ns.Experiment*ns.Epoch*ns.EpochChannel*ns.EpochParm &  channelRestrict & proj(tbl),'paradigm','signal','time','condition','align','ORDER BY channel');
            % Determine mean, ste and n averaging over the pv.average
            % dimension (anything that is left out of grouping)
            grouping = setdiff(["subject" "session_date" "starttime" "trial" "channel" "condition"  ],pv.average);
            funs  = {@(x) mean(cat(2,x{:}),2,"omitmissing")', ...  % Mean
                @(x) (std(cat(2,x{:}),0,2,"omitmissing")./sqrt(sum(~isnan(cat(2,x{:})),2,"omitmissing")))',...  % Standard error
                @(x) sum(~isnan(cat(2,x{:})),2,"omitmissing")'};  % Non-Nan N
            G = groupsummary(T, grouping, funs, 'signal');
            G= renamevars(G,["fun1_signal" "fun2_signal" "fun3_signal"],["mean" "ste" "n"]);
            % Combine with align/time/paradigm information. Note this
            % assumes these are constant across the group (picking
            % only the first or the unique here)
            P = groupsummary(T, grouping, @(x) x(1), "align");
            P = renamevars(P,"fun1_align","align");
            G =innerjoin(G,P);
            TM = groupsummary(T, grouping, @(x) unique(x)', ["time" "paradigm"]);
            TM  = renamevars(TM,["fun1_time" "fun1_paradigm"],["time" "paradigm"]);                        
            G =innerjoin(G,TM);
            G= sortrows(G,intersect(["subject" "session_date" "starttime" "condition" "channel" "trial"],G.Properties.VariableNames,'stable'));

            %% Figure
            tileCntr=0;
            nrTimeSeries = height(G);
            for i = 1:nrTimeSeries
                t =  G.time(i,:);
                t = linspace(t(1),t(2),t(3));
                align =G.align(i);
                m = G.mean(i,:);
                ste = G.ste(i,:);                
                n = mean(G.n(i,:));
                titlePV= setdiff(["paradigm" grouping],"condition");
                ttlStr = strjoin(string(G{i,titlePV}),"/");
                if i==1 || any(G{i,["subject" "session_date" "starttime"]} ~= G{i-1,["subject" "session_date" "starttime"]})
                    % New  subject, session or experiment in a new tile
                    if i>1
                        % Add the legend string to the existing tile
                        legend(h,legStr);
                    end
                    if mod(tileCntr,pv.tilesPerPage)==0
                        figure;
                    end
                    % Start a new tile with empty handles
                    nexttile;
                    tileCntr =tileCntr+1;
                    h = [];
                    legStr = string([]);
                    hold on
                end
                h = [h plot(t,m)];                %#ok<AGROW>
                patch([t flip(t)],[m+ste flip(m-ste)],h(end).Color,FaceAlpha= 0.5);
                legStr = [legStr G.condition(i)]; %#ok<AGROW>
                title (ttlStr + " (n=" + string(n) +")",'Interpreter','none');
                ylabel 'EP (\muV)'
                xlabel (sprintf('Time after %s.%s (ms)',align.plugin,align.event));
                % If delta is not empty, add the difference wave.
                if pv.delta ~="" && G.condition(i) ~=pv.delta
                    matchG = innerjoin(G,G(i,setdiff(grouping,"condition")));
                    reference = find(matchG.condition ==pv.delta);
                    if ~isempty(reference)
                        y = m - matchG.mean(reference,:);
                        ste = ste +matchG.ste(reference,:);
                        h = [h plot(t,y)];                     %#ok<AGROW>
                        patch([t flip(t)],[y+ste flip(y-ste)],h(end).Color,FaceAlpha= 0.5);
                        legStr = [legStr G.condition(i)+"-"+ pv.delta]; %#ok<AGROW>
                    end
                end
            end
            legend(h,legStr); % For the last tile

        end
    end
    methods (Access = protected)
        function makeTuples(tbl, key)

            %% Determine events to align to
            parmTpl = fetch(ns.EpochParm &key,'prep','art','align','window','channels');
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

            if isempty(parmTpl.channels)
                % Use all channels by default
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
            [signal,t,prepResults] = prep.preprocess(signal,seconds(t),parmTpl.prep,key);
            [nrSamples,nrTrials,nrChannels] = size(signal); %#ok<ASGLU>
            fprintf("\t Filtering is complete after %s\n",toc);
            %% --- Artifact/Outlier Rejection ---
            tic;
            fprintf("Artifact detection ...\n");
            samplingRate = 1./mode(diff(t));
            parmTpl.art.epoch_no = trials;
            pv =namedargs2cell(parmTpl.art);
            [artResult] = prep.artifactDetection(permute(signal,[2 3 1]),samplingRate,pv{:});
            % Remove epochs that were identified as having artifacts
            out = ismember(trials,artResult.all);
            signal(:,out,:) = [];
            trials(out) = [];
            condition(out) = [];
            startTime(out) = [];
            nrTrials =numel(trials);
            % Convert object to struct to allow storing in SQL database
            % with MYM
            warning('off','MATLAB:structOnObject')
            artResult = rmfield(struct(artResult),'categories');
            warning('on','MATLAB:structOnObject')
            fprintf("\t Artifact detection complete after %s\n",toc);

            %% --- Submit to the server ---
            tic;
            fprintf("Submitting to the server\n");
            epoch_tpl = mergestruct(key, ...
                struct(time = [t(1) t(end) numel(t)],...
                info = struct('prep',prepResults,'art',artResult)));
            % Insert to Epoch table
            chunkedInsert(tbl, epoch_tpl);

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
                emd
                fprintf("\t Submission is complete after %s.\n",toc);
            end
        end
    end
end

