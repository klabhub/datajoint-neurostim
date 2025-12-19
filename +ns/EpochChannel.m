%{
# Preprocessed and epoched data per channel and trial/epoch
-> ns.Epoch
channel : int       # Channel number
trial : int         # Trial number 
---
condition : varchar(128)  # Name of the condition to which this trial belongs
onset : float             # Time of the align event relative to trial start
signal : longblob         # C data for a single channel, single trial
%}
classdef EpochChannel < dj.Part & dj.DJInstance    
    properties (SetAccess = protected)
        master = ns.Epoch
        cache  (1,1) ns.epochCache;
    end

    properties (Dependent)
        channels                % Channels contributing to this Epoch table
    end

    methods      
        function ch = get.channels(self)
            ch = fetch(self, 'channel');
            ch = unique([ch(:).channel]');
        end
 
        function v = get.cache(tbl)
            % Avoid proj in the join  as this uses randomly assigned "AS"
            % names; as a result the cache will be refetched each time.               
            relvar = (ns.Experiment*ns.Epoch*ns.EpochChannel*ns.EpochParm) & proj(tbl);
            fill(tbl.cache,relvar); % Make sure the cache is filled and up to date.
            v = tbl.cache.T; % Return the table
        end
    end

    methods (Access = public)

      function plot (tbl,pv)
            % Plot evoked responses for all rows in the table.
            % Set the 'average' input to select which aspects to average
            % over. By default, trials and channels are averaged.
            %
            arguments
                tbl (1,1) ns.EpochChannel {mustHaveRows}
                pv.delta (1,1) string = ""              % Show the difference using this named condition as the reference
                pv.channel (:,1) double = []            % Select a subset of channels
                pv.trial (:,1) double = []            % Select a subset of trials
                pv.average (1,:) string {mustBeMember(pv.average,["starttime" "condition" "trial" "channel" "subject" "session_date" ""])} = ["trial" "channel"]  % Average over these dimensions
                pv.tilesPerPage (1,1) double = 6        % Select how many tiles per page.
                pv.linkAxes (1,1) logical = false        % Force the same xy axes on all tiles in a figure
                pv.raster (1,1) logical = false         % Set to true to show trials as rasters (removes "trial" from pv.average)
            end
            % Start a new tile when these values change:
            newTileEach = ["paradigm" "subject" "session_date" "starttime"];
            dimension = unique(tbl.cache.dimension);
            %% Fill the cache, then perform averaging per group            
            if pv.raster
                pv.average = setdiff(pv.average,"trial");
            end
            % Determine mean, ste and n averaging over the pv.average
            % dimension (anything that is left out of grouping)
            grouping = setdiff(["subject" "session_date" "starttime" "paradigm"  "condition" "trial" "channel"  ],pv.average,'stable');            
            G = compute(tbl.cache,@ns.EpochChannel.do_msten,["mean" "ste" "n"],channel=pv.channel,trial=pv.trial,grouping=grouping);
            
            if pv.raster
                % Concatenate the trials into a raster matrix in G.
                grouping = setdiff(grouping,"trial",'stable');
                P = groupsummary(G, grouping, @(x) x(1,:), ["align" "time"]);
                P = renamevars(P,["fun1_align" "fun1_time" ],["align" "time"]);
                G  = groupsummary(G,grouping,@(x) ({cat(1,x)}),["mean" "ste" "n"]);                
                G= renamevars(G,["fun1_mean" "fun1_ste" "fun1_n"],["mean" "ste" "n"]);                                       
                G =innerjoin(G,P);
                newTileEach = union(newTileEach,"condition");          
                G= sortrows(G,intersect([ "subject" "session_date" "starttime" "condition" "channel" "trial" "paradigm"],G.Properties.VariableNames,'stable'));
            end
                   
            
            %% Figure
            tileCntr=0;
            nrTimeSeries = height(G);
            newTileEach = intersect(newTileEach,G.Properties.VariableNames);
            legStr =string([]);
            for i = 1:nrTimeSeries
                t =  G.time(i,:);
                t = linspace(t(1),t(2),t(3));
                align =G.align(i);                
                if i==1 || (~isempty(newTileEach) && any(G{i,newTileEach} ~= G{i-1,newTileEach}))
                    % New  subject, session or experiment in a new tile
                    if i>1 &&  ~isempty(legStr)
                        % Add the legend string to the existing tile
                        legend(h,legStr);
                    end
                    if mod(tileCntr,pv.tilesPerPage)==0      
                        if i>1 && pv.linkAxes
                            linkaxes(gcf().Children().Children())
                        end
                        figure;
                    end
                    % Start a new tile with empty handles
                    nexttile;
                    tileCntr =tileCntr+1;
                    h = [];
                    legStr = string([]);
                    hold on
                end
                if pv.raster
                    % Show each condition in a separate tile
                    nrTrials= size(G.mean{i},1);
                    imagesc(t,1:nrTrials, G.mean{i})
                    n = mean(G.n{i},"all");
                    ylabel ("Trial")
                    titlePV= setdiff(["paradigm" grouping],"",'stable');
                    ttlStr = strjoin(string(G{i,titlePV}),"/");
                else
                    m = G.mean(i,:);
                    ste = G.ste(i,:);
                    n = mean(G.n(i,:));
                    h = [h plot(t,m)];                %#ok<AGROW>
                    p = patch([t flip(t)],[m+ste flip(m-ste)],h(end).Color,FaceAlpha= 0.5);                
                    p.EdgeColor = h(end).Color;
                    plot(xlim,[0 0],'k');
                    ylabel 'EP (\muV)'                    
                    legStr = [legStr dimension + "=" + G.condition(i)]; %#ok<AGROW>
                    titlePV= setdiff(["paradigm" grouping],"condition",'stable');
                    ttlStr = strjoin(string(G{i,titlePV}),"/");
                end
                title (ttlStr + " (n=" + string(n) +")",'Interpreter','none');                
                xlabel (sprintf('Time after %s.%s (s)',align.plugin,align.event));
                % If delta is not empty, add the difference wave.
                if pv.delta ~="" && G.condition(i) ~=pv.delta
                    matchVars = setdiff(grouping,"condition");
                    if isempty(matchVars)
                        matchG = G;
                    else    
                        matchG = innerjoin(G,G(i,matchVars));
                    end
                    reference = find(matchG.condition ==pv.delta);
                    if ~isempty(reference)
                        y = m - matchG.mean(reference,:);
                        ste = ste +matchG.ste(reference,:);
                        h = [h plot(t,y)];                     %#ok<AGROW>
                        p = patch([t flip(t)],[y+ste flip(y-ste)],h(end).Color,FaceAlpha= 0.5);
                        p.EdgeColor = h(end).Color;
                        legStr = [legStr G.condition(i)+"-"+ pv.delta]; %#ok<AGROW>
                    end
                end
            end
            if ~isempty(legStr); legend(h,legStr); end % For the last tile
            if pv.linkAxes
                linkaxes(gcf().Children().Children())
            end
        end


        function T = compute(self,fun, options,pv)
            % Compute derived measures from the EpochChannel table.
            % fun - fft, psd, pmtm, snr, msten
            % options - Struct passed to the compute function - different
            %           for each fun
            % channel  - Select a subset of channels
            % trial    - Select a subset of trials
            % grouping - Group analysis by these fields. Signals are
            %            averaged within group before applying the fun.
            arguments
                self (1,1) ns.EpochChannel
                fun  (1,1) string {mustBeMember(fun,["fft" "psd" "pmtm" "snr" "msten"])}
                options (1,:) struct = struct([]);  % Options for the fun
                pv.channel (:,1) double = []            % Select a subset of channels
                pv.trial (:,1) double = []            % Select a subset of trials
                pv.timeWindow (1,2) double = [-inf inf]  % Select a time window to operate on
                pv.grouping (1,:) string {mustBeMember(pv.grouping,["subject" "session_date" "starttime" "condition" "trial" "channel"])} = ["subject" "session_date" "starttime" "condition"]
            end
            eTbl = ns.Epoch & self;
            samplingRate = unique(round(eTbl.samplingRate));
            switch fun
                case "msten"
                    fun =@ns.EpochChannel.msten;
                    names = ["mean" "ste" "n"];
                case "fft"
                    assert(isempty(options),"fft does not take any options")
                    fun = @(x) ns.EpochChannel.do_fft(x,samplingRate);
                    names = ["fft" "phase" "frequency"];
                case "psd"
                   if isempty(options)
                        opts ={};
                   else
                        opts = namedargs2cell(options);
                   end
                    fun = @(x) ns.EpochChannel.do_psd(x,samplingRate,opts{:});
                    names = ["power" "frequency"];
                case "pmtm"                   
                    % opts.args = {4,7}
                    % opts.Tapers= 'Sine'
                    assert(isfield(options,"args"),"pmtm needs an .args field in the options struct, which will be passsed verbatim to pmtm.")
                    args =options.args;
                    opts = namedargs2cell(rmfield(options,'args')); % The rest as parm/value pairs
                    fun = @(x) ns.EpochChannel.do_pmtm(x,args{:},opts{:},samplingRate);
                    names = ["power" "frequency"];
                case "snr"
                    % First compute power using FFT
                    signalFun = @(x) ns.EpochChannel.do_fft(x,samplingRate);
                    signalNames = ["fft" "phase" "frequency"];
                    T = ns.EpochChannel.cacheCompute(self.cache,signalFun,signalNames ,channel =pv.channel,trial =pv.trial,timeWindow = pv.timeWindow,grouping=pv.grouping);

                    fq_res = mode(diff(T{1,"frequency"}));
                    n_bins = round(options.bin/fq_res);
                    n_bins_skip = round(options.bin_skip/fq_res);
                    mid_bin = n_bins + 1;
                    kernel = true(1, 2*n_bins + 1);
                    kernel(mid_bin-n_bins_skip : mid_bin+n_bins_skip) = false;
                    noiseFun = @(x) ns.EpochChannel.do_noise(x,kernel);
                    noiseNames = "noise";
                    N = ns.EpochChannel.cacheCompute(T,noiseFun,noiseNames,signal="fft",grouping=pv.grouping);
 
                    snr = cellfun(@(sig,nois) abs(sig)./nois,T.fft,N.noise,'UniformOutput',false);
                    T= addvars(T,snr);
                    T= removevars(T,signalNames);
                    return; % All done.
                otherwise 
                    error('Unknown function %s', fun);
            end
            % Send to cacheCompute
            T = ns.EpochChannel.cacheCompute(self.cache,fun,names,channel =pv.channel,trial =pv.trial,timeWindow = pv.timeWindow,grouping=pv.grouping);
        end
    end

   

    
   

end