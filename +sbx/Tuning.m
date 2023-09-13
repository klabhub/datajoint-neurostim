%{
# An ROI's direction tuning, fit with the sum of two von Mises functions.
-> sbx.Roi
-> ns.Experiment
-> sbx.TuningParms
-> ns.Condition
---
estimate   : blob # Tuning parameters [Preferred Peak AntiPeak Kappa Offset]
ci          : blob # Tuning parameters confidence intervals 
error       : blob # Error estimate
splithalves : float # Quality of the fit based on split-halves cross validation
r           : float # Correlation between fit and non-parametric estimate
residualz   : float # Z-score of mean residuals (0=good)
residualnoise:float # Ratio of stdev of residuals to the measurement noise (1 =good)
p           : float # Tuning p-value 
scale       : float  # Number by which the spike rate was scaled.
bootparms   : blob  # Boostrap estimates of the parameters
%}
%
% This table depends on sbx.Roi and sbx.Experiments, and the parameter
% settings for the Tuning table in the sbx.TuningParms table. 
% For instance:
%
% tuningParms  =struct('stimulus','gbr', ...            % Use the stimulus named gbr
%                        'independentVariable','orientation', ... % gbr.orientation has the directions/orientations
%                        'start',0.5,'stop',1.5,'step',1/15.5, 'interpolation','nearest', ...  % Time and bins within each trial to use
%                         'nrBoot',100, ...                              % Bootstrap resampling estimates
%                         'nrSplitHalves',20, ...                        % Splithalves validation (see poissyFit/splitHalves)
%                         'alpha',0.05, ...                          % Alpha level for CI
%                         'fPerSpike',500);                          % Scaling parameter "fluorescence per spike".
% Create a tuple, and tag this parameter set 'default' 
% tuningTpl = struct('tag','default','description','Fit to spikes fom 0.5 to 1.5 with two von Mises','parms',tuningParms);
% Insert into the sbx.TuningParms table
%    insert(sbx.TuningParms,tuningTpl)
% Then call 
% populate(sbx.Tuning)
%
% BK - Sept 2023
classdef Tuning <dj.Computed
    methods (Access=public)
        function plot(tbl,pv)
            arguments 
                tbl (1,1) sbx.Tuning
                pv.nrPerFigure (1,1) {mustBeNonnegative,mustBeInteger} = 20
                pv.nrBootToShow (1,1) {mustBeNonnegative,mustBeInteger} = 50
                pv.showRaw (1,1) logical =false
            end
            cntr =0;
            
            if pv.showRaw
             roi  = sbx.Roi & tbl;
             expt = ns.Experiment & tbl;
            end
            % Loop over the table
            for tpl = tbl.fetch('*')'
                cntr =cntr+1;
                % 1 figure for every nrPerFigure plots
                figure(ceil(cntr/pv.nrPerFigure));
                if mod(cntr-1,pv.nrPerFigure)+1==1
                    tiledlayout('flow')
                end
                nexttile;
                % Fetch parameters used to estmate the tuning curve
                parms = fetch1(sbx.TuningParms & tpl,'parms');
                
                % Generate estimated tuning curve
                tuningFunction = @poissyFit.logTwoVonMises;
                uStimulus = (0:1:360);
                
                if pv.nrBootToShow >0
                % Show Bootstrap estimates 
                bsParms = poissyFit.convertLogParms(tpl.bootparms,parms.step,"lin2log");                     
                hold on
                for i=1:min(pv.nrBootToShow,parms.nrBoot)
                    errorCurve =  exp(tuningFunction(uStimulus,bsParms(i,:)))'/parms.step;
                    plot(uStimulus,errorCurve,'LineWidth',0.5,'Color',0.6*ones(1,3),'LineStyle','-');
                end
                end
                
                % Prediction on top
                vonMisParms = poissyFit.convertLogParms(tpl.estimate',parms.step,"lin2log");              
                predictedTuningCurve = exp(tuningFunction(uStimulus,vonMisParms))'/parms.step;                
                plot(uStimulus,predictedTuningCurve,'LineWidth',2,'Color','k','LineStyle','-');

                % Show non parametric tuning cuvrev (mean/ste)
                if pv.showRaw
                     [~,spk] = get(roi &tpl,expt&tpl,trial = trials,modality = 'spikes',start=parms.start,stop=parms.stop,step=parms.step,interpolation = parms.interpolation);
                     spk = mean(spk,1,"omitnan");
                     direction = get(expt &tpl,parms.stimulus,'prm',parms.independentVariable,'atTrialTime',0);
                     [uDirection,~,ix] = unique(direction);
                     ix =repmat(ix',[size(spk,1) 1]);
                     tc = accumarray(ix(:),spk(:),[],@(x) mean(x,'all','omitnan'));
                     tcErr =accumarray(ix(:),spk(:),[],@(x) std(x,0,'all',"omitnan")./sqrt(sum(~isnan(x(:)))));                        
                     ploterr(uDirection, tc/parms.fPerSpike,tcErr/parms.fPerSpike,'ShadingAlpha',0.5,'LineWidth',2);
                end
                % Show parameters in the title
                 txt=  sprintf('Roi#%d - PD:%.0f [+/- %.0f], kappa: %.1f  [%.0f %.0f]\n Amp: %.2f [%.2f %.2f] AntiAmp: %.2f [%.2f %.2f], Offset %.2f  [%.2f %.2f] \n (p : %.3g)', ...
                        tpl.roi,...                      
                         tpl.estimate(1),tpl.error(1),...
                         tpl.estimate(4),tpl.ci(4,1),tpl.ci(4,2),...
                         tpl.estimate(2),tpl.ci(2,1),tpl.ci(2,2),...
                         tpl.estimate(3),tpl.ci(3,1),tpl.ci(3,2),...
                         tpl.estimate(5),tpl.ci(5,1),tpl.ci(5,2),...
                         tpl.p);
                 title(txt);
                 ylabel 'Spike Rate (spk/s)' 
                 drawnow;
            end
    
        end
    end
    methods (Access=protected)
        function makeTuples(tbl,key)
            if isempty(which('poissyFit'))
                error('This function uses the poissyFit class. Please install it first (https://github.com/klabhub/poissyFit)')
            end

            %% Retrive sources
            roi  = sbx.Roi & key;
            expt = ns.Experiment & key;
            parms = fetch1(sbx.TuningParms &key,'parms');
            trials = [fetch(ns.Condition*ns.ConditionTrial & key,'trial').trial];
              
                        
            if count(ns.Plugin & expt & struct('plugin_name',parms.stimulus))==0
                fprintf('This experiment does not have a %s  plugin (%d trials)\n. No tuning computed.',parms.stimulus, fetch1(expt,'trials'))
                return
            end
         
            
            
            % Get the data and the directions
            [~,spk] = get(roi,expt,trial = trials,modality = 'spikes',start=parms.start,stop=parms.stop,step=parms.step,interpolation = parms.interpolation);
            direction = get(expt ,parms.stimulus,'prm',parms.independentVariable,'atTrialTime',0);
            direction = direction(trials);
            % Remove trials with NaN
            out = any(isnan(spk),1);
            spk(:,out)= [];
            direction(out) =[];
            
            %% Fit
            % Setup the poissyFit function for spike rate estimation from
            % the deconvolved spike estimates; the "rate" will be
            % downscaled to match the maximum in the poissyFit class (100 spk/s).
            o = poissyFit(direction,spk,parms.step,@poissyFit.logTwoVonMises, ...
                            "fPerSpike",parms.fPerSpike, ...
                            "tau",eps, ...
                            "hasDerivatives",1, ...
                            "scaleToMax",false);
            o.spikeCountDistribution = 'POISSON';
            o.measurementNoise =estimateNoise(o);
            % If workers have been started, we'll use them
            if isempty(gcp('nocreate'))
                o.nrWorkers = 0;
            else
                o.nrWorkers = gcp('nocreate').NumWorkers;
            end                        
            o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
                'SpecifyObjectiveGradient',true, ...
                'display','none', ...
                'CheckGradients',false, ... % Set to true to check supplied gradients against finite differences
                'diagnostics','off');

            % Bootstrap the MLE of the tuning parms
            solve(o,parms.nrBoot);
            % Do split halves validation
            splithalves = splitHalves(o,parms.nrSplitHalves);
            
            % Translate parameters to more meaningful vars with CI that no
            % longer use the exp(parms) of poissyFit            
            %% Insert in table.
            tpl = mergestruct(key, ...
                struct('estimate',poissyFit.convertLogParms(o.parms,o.binWidth,"log2lin")', ...
                'ci',poissyFit.convertLogParms(o.parmsCI,o.binWidth,"log2lin")',...
                'error',poissyFit.convertLogParms(o.parmsError,o.binWidth,"log2lin")', ...
                'splithalves',splithalves, ...
                'r',o.gof, ...
                'residualz',o.residualZ, ...
                'residualnoise',o.residualNoise, ...
                'scale',o.scale, ...
                'p',o.p, ...
                'bootparms',poissyFit.convertLogParms(o.bootParms,o.binWidth,"log2lin")));
            insert(tbl,tpl);
        end
    end
end