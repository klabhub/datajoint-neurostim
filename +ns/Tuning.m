%{
# A tuning curve: signal averageds computed channel for each ns.Dimension
-> ns.CChannel
-> ns.TuningParm
-> ns.Dimension
---
peak      : float # peak value across conditions
preferred : float # condition/independent variable where the peak value occurs
mean      : float # Mean value over all conditions
min       : float # Minimum va;ie across all conditions.
tc        : blob  # Non parametric tuning curve (onve value per condition)
tcstd     : blob # Non parametric tuning curve standard deviation
tcx       : blob # Independent variable for the tc (e.g., orientation/direction)
panova     : float # ANOVA based on non parametric tuning curve.
baseline  : float # Mean baseline (spontaneous) value
baselinesd: float # Standard deviation of the baseline va;ie
nrtrials    : float # Number of trials used in estimates
nrtrialspercondition  : blob # Number of trials used per condition in estimates
fit=NULL    : blob # Parameters of a user-defined fit.  
%}
%
% The ns.TuningParms table defines how the tuning curves are computed.
% Each row in ns.TuningParms identifies a dimension (a row in ns.Dimension)
% for which the tuning is computed. In addition, the ns.TuningParms tuple
% contains a 'parms' struct with the following properties:
% (These are passed to ns.C/align)
% .start    The time in milliseconds at which to start averaging
% .stop     The time in milliseconds at which to stop averaging)
% .baselineStart The time in milliseconds at which to start averaging to determine the baseline
% .baselineStop The time in milliseconds at which to stop averaging to determine the baseline
% .interpolation How to average ('mean','median','sum','mode','min','max', or a function handle)
% .align  The name of an plugin that serves as t=0. If empty, times are
% defined relative to the first frame in each trial. Otherwise, start/stop
% are defined relative to the stimulus' startTime property.
%
% EXAMPLE
% To determine the mean response to a stimulus called 'gabor' in a
% window fomr 50 ms until 1050 ms after stimulus onset:
% parms.start  =50
% parms.stop = 1050;
% parms.baselineStart= 1100
% parms.baselineStop = 1200
% parns.interpolation = 'mean'
% parms.align = 'gabor'
%
% OPTIONAL  - Fit a parametric tuning curve
% .fun - A handle to a function with the prototype:
% out = fun(x,y,parms,npEstimate);
% out - a struct with results.  This struct will be stored in the table.
%       If this struct has a field .estimate  (the estimated parameter values) 
%       and .fun the function that takes the indepependent variable and the estimated 
%       parameter as input, then the plot funciton can show the results.
% This function will be called with
% x = The independent variable per trial
% y= The mean signal per trial
% parms  = The complete parms struct in ns.TuningParm (use this to
%           specify additional parameters for the fit).
% npEstimate - The nonparametric estimates (i.e. all values that will be
%               stored in the table, like .peak, .preferred,etc.). These
%               can, for instance, be used in your function to start the
%               fitting procedure from reasonable initial values.
% For an example, see ns.directionTuning.m
%
% BK - Nov 2023
classdef Tuning <dj.Computed
    properties (Dependent)
        keySource
    end

    methods

        function v = get.keySource(~)
            % Restricted to Dimensions listedin TuningParm
            dimTpl= [];
            for thisPrm = fetch(ns.TuningParm,'dimension')'
                restrict  =struct('dimension',thisPrm.dimension);
                thisTpl = fetch((ns.Dimension & restrict));
                if isempty(dimTpl)
                    dimTpl = thisTpl;
                else
                    dimTpl  = catstruct(1,dimTpl,thisTpl);
                end
            end
            % And then restrict the full table by the set of found tuples.
            v = (ns.CChannel*ns.TuningParm*proj(ns.Dimension &dimTpl)) ;
        end
    end

    methods (Access=public)
        function plot(tbl,pv)
            % Function to show a set of tuning functions
            arguments
                tbl (1,1) ns.Tuning
                pv.nrPerFigure (1,1) {mustBeNonnegative,mustBeInteger} = 20  % How many curves per figure
                pv.nrBootToShow (1,1) {mustBeNonnegative,mustBeInteger} = 50  % How many boostrap samples to show
                pv.linkAxes (1,1) logical =false                        % Link axes in the figure?
                pv.average (1,1) logical =false
            end
            cntr =0;
            % Loop over the table
            if pv.average
                nrChannels = 1;
            else
                nrChannels = count(tbl);
            end

            for roi =1:nrChannels
                if pv.average
                    tpl = fetch(tbl,'*');
                    x = tpl(1).tcx';
                    tc = cat(1,tpl.tc);
                    pk = cat(1,tpl.peak);
                    tc = tc./pk;
                    y = mean(tc,1)';
                    n = size(tc,1);
                    e = std(tc,0,1)'./sqrt(n);
                    base= cat(1,tpl.baseline)./pk;
                    b = mean(base);
                    be = std(base)./sqrt(n);
                else
                    tpl =fetch(tbl,'*',sprintf('LIMIT 1 OFFSET %d',roi-1));
                    x =tpl.tcx';
                    y =tpl.tc';
                    e = tpl.tcstd'./sqrt(numel(tpl.nrtrialspercondition));
                    b = tpl.baseline;
                    be = tpl.baselinesd/sqrt(tpl.nrtrials);
                end
                cntr =cntr+1;
                % 1 figure for every nrPerFigure plots
                figure(ceil(cntr/pv.nrPerFigure));
                if mod(cntr-1,pv.nrPerFigure)+1==1
                    T= tiledlayout('flow');
                end
                nexttile;
                hold on


                % The .fit field is optional for user defined parametric
                % fits and can have any number of fields, this is
                % specifically for the struct output of
                % ns.directionTuning.m although other fits could follow the
                % same format.
                if ~pv.average && isfield(tpl,'fit') 
                    % Show estimated fitted tuning curve
                    if pv.nrBootToShow >0 && isfield(tpl.fit,'bootparms')
                        % Show Bootstrap estimates
                        for i=1:min(pv.nrBootToShow,tpl.parms.nrBoot)
                            curve =  feval(tpl.fit.fun,tpl.tcx,tpl.fit.bootparms(i,:))';
                            plot(tpl.tcx,curve,'LineWidth',0.5,'Color',0.6*ones(1,3),'LineStyle','-');
                        end
                    end
                    % Prediction on top
                    if isfield(tpl.fit,'estimate')
                        predictedTuningCurve = feval(tpl.fit.fun,tpl.tcx,tpl.fit.estimate)';
                        plot(tpl.tcx,predictedTuningCurve,'LineWidth',2,'Color','k','LineStyle','-');
                    end
                end

                % Show non parametric tuning cuvrev (mean/ste)
                % Get the data and the directions
                ploterr(x,y,e,'ShadingAlpha',0.5,'LineWidth',2,'Color','r');
                ploterr(xlim',b*[1 1]',be*[1 1]','ShadingAlpha',0.5,'Color','k');
                % Show parameters in the title
                if pv.average
                else
                    txt=  sprintf('%s (%d trials)- Roi#%d \n', ...
                        tpl.tuningtag,...
                        tpl.nrtrials,...
                        tpl.channel);
                    title(txt,'Interpreter','None');
                end

                ylabel 'Spike Rate (spk/s)'
                drawnow;
                if pv.linkAxes
                    linkaxes(T.Children)
                end
            end

        end
    end
    methods (Access=protected)

        function makeTuples(tbl,key)
            tic
            %% Retrieve sources
            parms = fetch1(ns.TuningParm &key,'parms');
            conditions = fetch(ns.DimensionCondition & key,'trials','value');
            assert(numel(conditions)>0,"No conditions in %s group",key.dimension)
            %% Find indepdent variable values from the name assigned to the condition (which is
            % pluginNane:parameterName:Value. Assuming these values are
            % double..

            % Get the spikes for all trials (some may not be part of this
            % condition group)
            [spk ,~]= align(ns.C & key,channel =key.channel, start=parms.start,stop=parms.stop,step=(parms.stop-parms.start),interpolation = parms.interpolation);
            spk =spk{1,:};
            xValue= [conditions.value];
            if iscell(xValue)
                xValue = [xValue{:}];
            end
            % Sort by x
            [xValue,ix] = sort(xValue);
            conditions =conditions(ix);
            trialsPerCondition = {conditions.trials};

            % Average
            tuningCurve = cellfun(@(x) mean(spk(x),"omitnan"),trialsPerCondition);
            tuningCurveSd = cellfun(@(x) std(spk(x),0,"omitnan"),trialsPerCondition);
            nrTrialsPerCondition = cellfun(@numel,trialsPerCondition);
            [peak,ix] = max(tuningCurve);
            preferred = xValue(ix);
            m = mean(tuningCurve);
            mi = min(tuningCurve);

            % Map trials to conditions to do anova.
            conditionIx = nan(1,size(spk,2));
            for i=1:numel(conditions)
                conditionIx(conditions(i).trials)=i;
            end
            stayTrials = ~isnan(conditionIx);
            conditionIx(~stayTrials) = [];
            pAnova = anovan(spk(stayTrials),conditionIx','display','off');
            if isnan(pAnova);pAnova=1;end

            % Extract baseline (use only trials that are alos used for the
            % tuning).
            [spkBaseline,~] = align(ns.C&key,channel=key.channel,start=parms.startBaseline,stop=parms.stopBaseline,step=(parms.stopBaseline-parms.startBaseline),interpolation = parms.interpolation);
            spkBaseline = spkBaseline{1,:};
            baseline = mean(spkBaseline(stayTrials),"all","omitnan");
            baselineStd = std(spkBaseline(stayTrials),0,"all","omitnan");

            % Combine results into a struct
            estimate = struct('peak',peak,...
                'mean',m,...
                'preferred',preferred,...
                'tc',tuningCurve,...
                'tcstd',tuningCurveSd,...
                'tcx',xValue,...
                'min',mi,...
                'panova',pAnova, ...
                'nrtrialspercondition',nrTrialsPerCondition, ...
                'nrtrials',sum(nrTrialsPerCondition),...
                'baseline',baseline, ...
                'baselinesd',baselineStd);

            % If fun is specified; do a parametric fit using the specified function
            if isfield(parms,'fun')
                x = xValue(conditionIx);
                y = spk(stayTrials);
                out = isnan(y); % Remove trials with NaN estimates (happens for first trial sometimes)
                x(out)=[];
                y(out)=[];
                estimate.fit = feval(parms.fun,x,y,parms,estimate);
            end


            %% Insert in table.
            tpl = mergestruct(key, estimate);
            insert(tbl,tpl);

        end
    end
end