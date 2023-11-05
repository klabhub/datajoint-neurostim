%{
# An ROI's nonparametreic tuning: rates per condition 
-> sbx.Roi
-> ns.Experiment
-> sbx.TuningParms
---
peak      : float # peak firing rate
preferred : float # prefferred condition/independent variable.
mean      : float # Mean firing rate over all conditions
min       : float # Minimum response across all conditions.
tc        : blob  # Non parametric tuning curve (onve value per condition)
tcstd     : blob # non parametric tuning curve standard deviation
tcx       : blob # Independent variable for the nptc (e.g., orientation/direction)
panova     : float # ANOVA based on non parametric tuning curve.
baseline  : float # Mean baseline (spontaneous) activity
baselinesd: float # Standard deviation of the baseline activity
nrtrials    : float # Number of trials used in estimates
nrtrialspercondition  : blob # Number of trials used in estimates
%}
%
% This table depends on sbx.Roi and sbx.Experiments, and the parameter
% settings in the sbx.NPTuningParms table.
%
% BK - Nov 2023
classdef NPTuning <dj.Computed
    properties (Dependent)
        keySource
    end
    methods 
        function v= get.keySource(~)
            % Use only those experiments that have an entry in the Trialmap
            % (as that means that the sbx data were processed and the time course of 
            % the ROI is known in the experiment). Basically this makes
            % sure we ignore experiments where the imaging data were
            % missing.            
            v = (ns.Experiment & sbx.PreprocessedTrialmap)*sbx.Roi*sbx.TuningParms;            
        end
    end


    methods (Access=public)
        function plot(tbl,pv)
            % Function to show a set of tuning functions
            arguments
                tbl (1,1) sbx.NPTuning
                pv.nrPerFigure (1,1) {mustBeNonnegative,mustBeInteger} = 20  % How many curves per figure
                pv.nrBootToShow (1,1) {mustBeNonnegative,mustBeInteger} = 50  % How many boostrap samples to show
                pv.showNP (1,1) logical =false                          % Show Non-Parametric estimates?
                pv.linkAxes (1,1) logical =false                        % Link axes in the figure?
                pv.average (1,1) logical =false
            end
            cntr =0;
            % Loop over the table
            if pv.average                 
                nrRoi = 1;
            else
                nrRoi = count(tbl);
            end

            for roi =1:nrRoi
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
                    numel(tpl.nrtrials),...
                    tpl.roi);
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
            roi  = sbx.Roi & key;
            expt = ns.Experiment & key;
            parms = fetch1(sbx.TuningParms &key,'parms');
            conditions = fetch(ns.Condition & key & ['condition_group=''' parms.cgroup ''''],'trials');
            assert(numel(conditions)>0,"No conditions in %s group",parms.cgroup)
            %% Find indepdent variable balues from the name assigned to the condition (which is 
            % pluginNane:parameterName:Value. Assuming these values are
            % double..
            xValue= str2double(extractAfter({conditions.name},[parms.independentVariable ':']));
            % Sort by x
            [xValue,ix] = sort(xValue);
            conditions =conditions(ix); 
            % Get the spikes for all trials (some may not be part of this
            % condition group)
            [~,spk] = get(roi,expt,modality = 'spikes',start=parms.start/1000,stop=parms.stop/1000,step=(parms.stop-parms.start)/1000,interpolation = parms.interpolation);
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
            [~,spkBaseline] = get(roi,expt,modality = 'spikes',start=parms.startBaseline/1000,stop=parms.stopBaseline/1000,step=(parms.stopBaseline-parms.startBaseline)/1000,interpolation = parms.interpolation);
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