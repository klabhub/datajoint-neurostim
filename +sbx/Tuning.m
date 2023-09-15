%{
# An ROI's direction tuning, fit with the sum of two humps (von Mises functions.
-> sbx.Roi
-> ns.Experiment
-> sbx.TuningParms
-> ns.Condition
---
amplitude   : float # Response at the preferred direction. 
preferred   : float # Preferred direction in degrees
width       : float # Widht of the tuning at the preferred
antiwidth   : float # Widht of the tuning at the anti preferred
antiamplitude:float # Response at the anti-preferred direction.
offset      : float # Offset response.
ci          : blob # Tuning parameters confidence intervals  [Offset Preferred Width AntiWidth Peak AntiPeak]
error       : blob # Error estimate [Offset Preferred Width AntiWidth Peak AntiPeak]
splithalves : float # Quality of the fit based on split-halves cross validation
r           : float # Correlation between fit and non-parametric estimate
residualz   : float # Z-score of mean residuals (0=good)
bootparms   : blob  # Boostrap estimates of the parameters
nppeak      : float # Non parameteric peak firing rate
nppreferred : float # non parametric prefferred direction.
npmean      : float # Mean firing rate over all directions
npmin       : float # Minimum response.
nptc        : blob  # Non parametric tuning curve
nptcerr     : blob # non parametric tuning curve standard errror
nptcx       : blob # Independent variable for the nptc (orientation/direction)
npanova     : float # ANOVA based on non parametric tuning curve.
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
%                         'alpha',0.05);                          % Alpha level for CI
%
% Create a tuple, and tag this parameter set 'default'
% tuningTpl = struct('tag','default','description','Fit to spikes fom 0.5 to 1.5 with two humps','parms',tuningParms);
% Insert into the sbx.TuningParms table
%    insert(sbx.TuningParms,tuningTpl)
% Then call
% populate(sbx.Tuning)
%
%
% Although this table is specific for direction tuning, it can easily be
% extended to fit other parametric functions too.
%
% BK - Sept 2023
classdef Tuning <dj.Computed
    methods (Access=public)
        function plot(tbl,pv)
            arguments
                tbl (1,1) sbx.Tuning
                pv.nrPerFigure (1,1) {mustBeNonnegative,mustBeInteger} = 20
                pv.nrBootToShow (1,1) {mustBeNonnegative,mustBeInteger} = 50
                pv.showNP (1,1) logical =false
                pv.linkAxes (1,1) logical =false
            end
            cntr =0;


            % Loop over the table
            for tpl = tbl.fetch('*')'
                cntr =cntr+1;
                % 1 figure for every nrPerFigure plots
                figure(ceil(cntr/pv.nrPerFigure));
                if mod(cntr-1,pv.nrPerFigure)+1==1
                    T= tiledlayout('flow');
                end
                nexttile;
                % Fetch parameters used to estmate the tuning curve
                parms = fetch1(sbx.TuningParms & tpl,'parms');
                cond = fetch(ns.Condition*ns.ConditionTrial & tpl,'*');
                conditionName = unique({cond.name});
                conditionName= conditionName{1};
                trials = [cond.trial];

                % Generate estimated tuning curve
                uDirection = (0:1:330);                
                if pv.nrBootToShow >0
                    % Show Bootstrap estimates
                    hold on
                    for i=1:min(pv.nrBootToShow,parms.nrBoot)
                        errorCurve =  sbx.Tuning.twoHumps(uDirection,tpl.bootparms(i,:))';
                        plot(uDirection,errorCurve,'LineWidth',0.5,'Color',0.6*ones(1,3),'LineStyle','-');
                    end
                end

                % Prediction on top
                estimate = [tpl.offset tpl.preferred tpl.width tpl.antiwidth tpl.amplitude tpl.antiamplitude];
                predictedTuningCurve = sbx.Tuning.twoHumps(uDirection,estimate)';
                plot(uDirection,predictedTuningCurve,'LineWidth',2,'Color','k','LineStyle','-');

                % Show non parametric tuning cuvrev (mean/ste)
                if pv.showNP
                    % Get the data and the directions
                    ploterr(tpl.nptcx, tpl.nptc,tpl.nptcerr/sqrt(numel(trials)),'ShadingAlpha',0.5,'LineWidth',2,'Color','r');
                end
                % Show parameters in the title
                txt=  sprintf('%s (%d trials)- Roi#%d \n PD:%.0f [+/- %.0f], Amp: %.2g [%.2g %.2g], AntiAmp: %.2g [%.2g %.2g], Width: %.2g  [%.2g %.2g]  ,AntiWidth: %.2g  [%.2g %.2g]  , Offset %.2g  [%.2g %.2g] \n gof : %.3f, splithalves: %.2f', ...
                    conditionName, ...
                    numel(trials),...
                    tpl.roi,...
                    tpl.preferred,tpl.error(2),...
                    tpl.amplitude,tpl.ci(5,1),tpl.ci(5,2),...
                    tpl.antiamplitude,tpl.ci(6,1),tpl.ci(6,2),...
                    tpl.width,tpl.ci(3,1),tpl.ci(3,2),...
                    tpl.antiwidth,tpl.ci(4,1),tpl.ci(4,2),...
                    tpl.offset,tpl.ci(1,1),tpl.ci(1,2),...
                    tpl.r, ...
                    tpl.splithalves);
                title(txt);
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
            %% Retrive sources
            roi  = sbx.Roi & key;
            expt = ns.Experiment & key;
            parms = fetch1(sbx.TuningParms &key,'parms');
            trials = [fetch(ns.Condition*ns.ConditionTrial & key,'trial').trial];
            if count(ns.Plugin & expt & struct('plugin_name',parms.stimulus))==0
                fprintf('This experiment does not have a %s  plugin (%d trials)\n. No tuning computed.',parms.stimulus, fetch1(expt,'trials'))
                return
            end
tic
            % Get the data and the directions
            [~,spk] = get(roi,expt,trial = trials,modality = 'spikes',start=parms.start,stop=parms.stop,step=parms.step,interpolation = parms.interpolation);
            direction = get(expt ,parms.stimulus,'prm',parms.independentVariable,'atTrialTime',0);
            direction = direction(trials);
            % Remove trials with NaN
            out = any(isnan(spk),1);
            spk(:,out)= [];
            direction(out) =[];
            % If multiple bins per trial are used, we need to repmat the
            % direction to match those bins  (in that case each bin is
            % considered an independent observation to be fitted). If you
            % just want to fit the mean rate in a trial, set step to
            % stop-start which will generate a single bin per trial.
            direction = repmat(direction,[size(spk,1) 1]);

            spk = spk(:); % single column
            direction = direction(:); % Column


            %% Determine the nonparametric tuning curve
            [uDirection,~,stimulusIx] = unique(direction);
            tuningCurve = accumarray(stimulusIx,spk,[],@(x) mean(x,"omitnan"));
            tuningCurveErr = accumarray(stimulusIx,spk,[],@(x) std(x,0,1,"omitnan"));
            [npAmplitude,ix] = max(tuningCurve);
            npPreferred = uDirection(ix);
            npMean = mean(tuningCurve);
            npAntiAmplitude  = tuningCurve(uDirection==npPreferred+180 | uDirection==npPreferred-180);
            npMin = min(tuningCurve);
            npAnova = anovan(spk,direction,'display','off');


            %% Fit a parametric (twoHumps) function
            bootfun = @(x,y) sbx.Tuning.solve(x,y,npMin,npPreferred,npAmplitude,npAntiAmplitude);
            bootOpts = statset;
            bootOpts.UseParallel = true;
            [bootEstimates]=bootstrp(parms.nrBoot,bootfun,direction,spk,'Options',bootOpts);
            % Determine preferred axis and circular standard deviation
            % Multiply by two in case the preferred swaps
            % between bootstrap sets (i.e. an orientation tuned
            % neuron).
            z = exp(2*pi/180*1i*bootEstimates(:,2));
            R = mean(z);

            preferredOrAntiPreferred = angle(R)/2; % Back to the original 360 space
            % Because we always call the direction with the highest amplitude the
            % preferred direction (see "SOLVE" above), we just have to count which occurs
            % most frequently : preferredAxis or preferredAxis+180.
            z = exp(pi/180*1i*bootEstimates(:,2));
            inner = [real(z), imag(z)]*[cos(preferredOrAntiPreferred); sin(preferredOrAntiPreferred)];
            if mean(inner>0)>0.5
                isPreferred = inner>0;
            else
                isPreferred = inner<0;
            end

            bootEstimates(~isPreferred,[2 4 5]) = bootEstimates(~isPreferred,[2 5 4]) + [180 0 0];
            bootEstimates(:,2)= mod(bootEstimates(:,2),360);% Convenience

            estimate = mean(bootEstimates,1,'omitnan');
            error = std(bootEstimates,0,1,"omitnan");
            ci= prctile(bootEstimates,[parms.alpha/2 1-parms.alpha/2]);

            % Correct circular vriables.
            z = exp(pi/180*1i*bootEstimates(:,2));
            R = mean(z);
            estimate(2) = 180/pi*angle(R);
            error(2) = 180/pi*sqrt(-2*log(abs(R))); % Circular standard deviation
            % Use circular standard deviation as the
            % confidence limits.
            ci(:,2) = mod(estimate(2) + error(2)*[-1 1]',360);

            fittedCurve = sbx.Tuning.twoHumps(uDirection,estimate);
            gof = corr(tuningCurve,fittedCurve,'type','Spearman');

            residuals = sbx.Tuning.twoHumps(direction,estimate)- spk;
            meanR = mean(residuals,"all","omitnan");
            stdR =  std(residuals,0,"all","omitnan");
            residualZ = abs(meanR./stdR);

            %% SplitHalves
            r= nan(1,parms.nrSplitHalves);
            for i=1:parms.nrSplitHalves
                [oneHalfTrials,otherHalfTrials] =resampleTrials(stimulusIx,false,0.5) ;
                % First half
                firstHalfEstimate = sbx.Tuning.solve(direction(oneHalfTrials),spk(oneHalfTrials),npMin,npPreferred,npAmplitude,npAntiAmplitude);
                otherHalfEstimate = sbx.Tuning.solve(direction(otherHalfTrials),spk(otherHalfTrials),npMin,npPreferred,npAmplitude,npAntiAmplitude);
                tc1 = sbx.Tuning.twoHumps(uDirection,firstHalfEstimate);
                tc2 = sbx.Tuning.twoHumps(uDirection,otherHalfEstimate);
                r(i) = corr(tc1,tc2);
            end
            splithalves  = mean(r);


            %% Insert in table.
            tpl = mergestruct(key, ...
                struct('offset',estimate(1),...
                'preferred',estimate(2),...
                'width',estimate(3),...
                'antiwidth',estimate(4),...
                'amplitude',estimate(5), ...
                'antiamplitude',estimate(6),...
                'ci',ci',...
                'error',error', ...
                'splithalves',splithalves, ...
                'r',gof,...
                'residualz',residualZ,...
                'bootparms',bootEstimates,...
                'nppeak',npAmplitude,...
                'npmean',npMean,...
                'nppreferred',npPreferred,...
                'nptc',tuningCurve,...
                'nptcerr',tuningCurveErr,...
                'nptcx',uDirection,...
                'npmin',npMin,...
                'npanova',npAnova));
            insert(tbl,tpl);
            toc
        end
    end

    methods (Static)
        function [y] = twoHumps(x,parms)
            % Use this to fit a direction selective tuning function with
            % one bump at the preferred, and a (potentially smaller) bump
            % at the anti-preferred (=preferred +180).
            %
            % I tried this with the sum of two Von Mises functions:
            % y = amp1*exp(kappa1*cos((x - preferred))) +amp2*exp(kappa2*cos((x - preferred-180)))+offset;
            % but did not work well (solutions were oddly always pretty
            % poor). Currently the function uses two wrapped/circular
            % gaussians instead, which seems to work better.
            %
            % INPUT
            % x - Stimulus angles in degrees.
            % parms - parameter vector - [offset, preferred , width1, width2, amp1 ,amp2]
            %
            % OUTPUT
            %  y - Value per x
            offset = parms(1); preferred = parms(2); sigma1 = parms(3); sigma2 = parms(4); amp1=parms(5); amp2=parms(6);

            % Von Mises attempt
            %deg2rad =pi/180;
            %term1 = amp1*exp(sigma1*cos(deg2rad*(x-preferred)));
            %term2 = amp2*exp(sigma2*cos(deg2rad*(x-preferred-180)));
            %y = term1 + term2 + offset;

            y = offset;
            for k=-4:4
                y = y + amp1*exp(-(1/sigma1^2)*(x-preferred-360*k).^2);
                y = y + amp2*exp(-(1/sigma2^2)*(x-preferred-180-360*k).^2);
            end
        end

        function estimate = solve(direction,spk,npMin,npPreferred,npAmplitude,npAntiAmplitude)
            % Use the optimization toolbox to do a nonlinear least squares
            % fit of the spiking data to the twoHumps function

            % Setup the variables and bounds
            offset = optimvar('offset',1,'LowerBound',0,'UpperBound',2*npMin);
            preferred = optimvar('preferred',1,'LowerBound',npPreferred-45,'UpperBound',npPreferred+45);
            sigma = optimvar('sigma',1,'LowerBound',0);
            amplitude = optimvar('amplitude',1,'LowerBound',0,'UpperBound',2*npAmplitude);
            antiAmplitude = optimvar('antiAmplitude',1,'LowerBound',0,'UpperBound',2*npAntiAmplitude);
            antiSigma = optimvar('antiSigma',1,'LowerBound',0); % kappa= 30 is a 25 degree FWHM.

            % Least square
            objective = sum((spk  - sbx.Tuning.twoHumps(direction,[offset preferred sigma antiSigma amplitude antiAmplitude ])).^2);
            problem = optimproblem("Objective",objective);

            % Initial values
            initialValue.offset = 0;
            initialValue.sigma= 100;            % Start wide
            initialValue.antiSigma = 100;
            initialValue.amplitude= npAmplitude;  % Start with nonparameteric estiamtes
            initialValue.antiAmplitude = npAntiAmplitude;
            initialValue.preferred = npPreferred;

            % Set options
            opts =optimset;
            opts.Display = 'off';
            opts.Algorithm = 'trust-region-reflective';

            % Solve and return as vector
            solution =solve(problem,initialValue,'Options',opts);
            estimate  = [solution.offset mod(solution.preferred,360) solution.sigma  solution.antiSigma solution.amplitude solution.antiAmplitude];
        end
    end

end