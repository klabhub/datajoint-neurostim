%{
# A tuning curve: rates per condition 
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
%}
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
                pv.showNP (1,1) logical =false                          % Show Non-Parametric estimates?
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
                     y = mean(cat(1,tpl.tc),1)';
                     e = mean(cat(1,tpl.tcstd)./sqrt(numel([tpl.nrtrialspercondition])),1)';
                     b = mean(cat(1,tpl.baseline));
                     be = mean(cat(1,[tpl.baseline]./sqrt([tpl.nrtrials])));
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
            conditions = fetch(ns.DimensionCondition & key,'trials');
            assert(numel(conditions)>0,"No conditions in %s group",key.dimension)
            %% Find indepdent variable values from the name assigned to the condition (which is 
            % pluginNane:parameterName:Value. Assuming these values are
            % double..
            xValue= str2double(extractAfter({conditions.name},[parms.independentVariable ':']));
            % Sort by x
            [xValue,ix] = sort(xValue);
            conditions =conditions(ix); 
            % Get the spikes for all trials (some may not be part of this
            % condition group)
            [spk ,~]= align(ns.C & key,channel =key.channel, start=parms.start,stop=parms.stop,step=(parms.stop-parms.start),interpolation = parms.interpolation);
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

            % Extract baseline (use only trials that are alos used for the
            % tuning).
            [spkBaseline,~] = align(ns.C&key,channel=key.channel,start=parms.startBaseline,stop=parms.stopBaseline,step=(parms.stopBaseline-parms.startBaseline),interpolation = parms.interpolation);
            baseline = mean(spkBaseline(stayTrials),"all","omitnan");
            baselineStd = std(spkBaseline(stayTrials),0,"all","omitnan");


         
            %% Insert in table.
            tpl = mergestruct(key, ...
                struct('peak',peak,...
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
                'baselinesd',baselineStd));
            insert(tbl,tpl);

   
            toc

         
        end
    end
end