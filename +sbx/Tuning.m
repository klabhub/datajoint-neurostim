%{
# An ROI's direction tuning, fit with the sum of two von Mises functions.
-> sbx.Roi
-> ns.Experiment
---
parms       : blob # Tuning parameters [Preferred Peak AntiPeak Kappa Offset]
ci          : blob # Tuning parameters confidence intervals 
error       : blob # IQR of the parameters
splithalves : float # Quality of the fit based on split-halves cross validation
r           : float # Correlation between fit and non-parametric estimate
residualz   : float # Z-score of mean residuals (0=good)
residualnoise:float # Ratio of stdev of residuals to the measurement noise (1 =good)
stepsize    : float # Binwidth used
p           : float # Tuning p-value 
bootstrap   : integer # Number of bootstrap sets used
scale       : float  # Number by which the spike rate was scaled.
alpha       : float # Alpha level for CI
%}
classdef Tuning <dj.Computed

    methods (Access=protected)
        function makeTuples(tbl,key)
            if isempty(which('poissyFit'))
                error('This function uses the poissyFit class. Please install it first (https://github.com/klabhub/poissyFit)')
            end

            %% Setup and retrieve ROI data.
            roi  = sbx.Roi & key;
            expt = ns.Experiment & key;

            persistent cntr
            if isempty(cntr)
                cntr= 1;
            else
                cntr= cntr+1;
            end
            
            % Hardcoded, but stored in table.
            stepSize = 1/15.5;  % Always using the full sampling rate.
            nrBoot = 100;                                           
            show = false; % For debugging, set to true
            alpha = 0.05;

            % If workers have been started, we'll use them
            if isempty(gcp('nocreate'))
                nrWorkers = 0;
            else
                nrWorkers = gcp('nocreate').NumWorkers;
            end            
            % Get the data and the directions
            [~,spk] = get(roi,expt,modality = 'spikes',start=0.5,stop=1.5,step=stepSize,interpolation ='nearest');
            direction = get(expt ,'gbr','prm','orientation','atTrialTime',0);
            % Rmeove trials with NaN
            out = any(isnan(spk),1);
            spk(:,out)= [];
            direction(out) =[];
            
            %% Fit
            % Setup the poissyFit function for spike rate estimation from
            % the deconvolved spike estimates; the "rate" will be
            % downscaled to match the maximum in the poissyFit class (100 spk/s).
            o = poissyFit(direction,spk,stepSize,@poissyFit.logTwoVonMises, ...
                            "fPerSpike",1, ...
                            "tau",0, ...
                            "hasDerivatives",1, ...
                            "scaleToMax",true);
            o.spikeCountDistribution = 'POISSON';
            o.measurementNoise =estimateNoise(o);
            o.nrWorkers = nrWorkers;
            o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
                'SpecifyObjectiveGradient',true, ...
                'display','none', ...
                'CheckGradients',false, ... % Set to true to check supplied gradients against finite differences
                'diagnostics','off');

            % Bootstrap the MLE of the tuning parms
            solve(o,nrBoot);
            % Do split halves validation
            splithalves = splitHalves(o,nrBoot);
            
            % Translate parameters to more meaningful vars with CI that no
            % longer use the exp(parms) of poissyFit
            [~,prefPeak,antiPeak,kappa,offset] = poissyFit.twoVonMisesParms(o.bootParms,o.binWidth);
            bsParms = [prefPeak antiPeak kappa offset];
            ci    = prctile(bsParms,[100*alpha/2  100*(1-alpha/2)]);
            err = iqr(bsParms);
            parms  = median(bsParms);            
            parms = [mod(o.parms(2),360) parms]; % Add preferred 
            ci =  [o.parmsCI(:,2) ci];
            err = [o.parmsError(:,2) err];

            %% Visualizee
            % If debugging/visualizeing , show the tuning curves and
            % estimated parameters
            if show
                figure(ceil(cntr/20));
                if mod(cntr-1,20)+1==1
                    tiledlayout('flow')
                end
                nexttile;
                plot(o,showErrorbars= false, showBootstrapSets= true,equalAxes=true);                
                txt=  sprintf('Roi#%d - PD:%.0f [%.0f %.0f], kappa: %.1f  [%.0f %.0f]\n Amp: %.2f [%.2f %.2f] AntiAmp: %.2f [%.2f %.2f], \n (p : %.3g)',key.roi,...
                        parms(1),ci(1,1),ci(2,1),...
                        parms(4),ci(1,4),ci(2,4),...
                        parms(2),ci(1,2),ci(2,2),...
                        parms(3),ci(1,3),ci(2,3),...
                        o.p);
                title(txt);
                yyaxis left
                ylabel 'Spike Rate (spk/s)' % Avoid confusion
                legend('Data','Fit')
                drawnow;
            end
    
            %% Insert in table.
            tpl = mergestruct(key, ...
                struct('parms',parms', ...
                'ci',ci',...
                'error',err', ...
                'splithalves',splithalves, ...
                'r',o.gof, ...
                'residualz',o.residualZ, ...
                'residualnoise',o.residualNoise, ...
                'bootstrap',nrBoot, ...
                'stepsize',stepSize, ...
                'scale',o.scale, ...
                'p',o.p,...
                'alpha',alpha));
            insert(tbl,tpl);
        end
    end
end