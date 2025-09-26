%{ 
# Calibrated MLSpike parameters for a session.
-> sbx.Preprocessed
-> sbx.SpikesParm
---
quality = NULL: float             # Correlation between reconstructed F and F
a       = NULL: float           # Calibrated F per spike
sigma   = NULL: float               # Calibrated noise level    
tau     = NULL: float                 # Calibrated decay parameter
failed  : smallint           # Number of ROIs where calibration failed.
%}
classdef Mlspikecalibration < dj.Computed
    properties (Dependent)
        keySource
    end

    methods 
        function v = get.keySource(tbl)
            v = sbx.Preprocessed * (sbx.SpikesParm & 'calibration IS NOT NULL ');
        end
    end
    methods (Access=public)
         function G = plot(tbl,pv)
            % Show meta data on MLspike deconvolution. Autocal should be
            % successful on most rois and the quality of the reconstruction
            % (quality) should be high (r>0.8). Deviations from this
            % suggest that some parameters need to be tuned. The red lines
            % show the default values that will be used when autocal fails.
            % Ideally those would be near the center of the autocal
            % distributions. (Or one could be conservative and not analyze
            % any ROIs for which autocal fails). For instance with a query
            % like this :
            % ns.CChannel & proj(sbx.Spikes & 'autocal=1' &  'quality>0.8','roi->channel')
            arguments
                tbl (1,1) sbx.Spikes {mustHaveRows(tbl)}
                pv =1
            end

            for tpl = fetch(sbx.SpikesParm & tbl,'*')'

                T = fetchtable(tbl &tpl ,'*');
                figByName(tpl.stag)
                clf

                layout = tiledlayout('flow');
                nexttile
                % Quality (R2) for those where autocal worked /not worked
                bins = linspace(0,1,20);
                histogram(T.quality(T.autocal==1),bins);
                hold on
                histogram(T.quality(T.autocal==0),bins);
                title(sprintf('Quality (r = %.2f), autocal success: %.0f%%',mean(T.quality),100*mean(T.autocal)));
                xlabel 'r'
                ylabel '#roi'


                nexttile
                % Noise level  Sigma for those where autocal worked /not worked
                bins = linspace(0,1,20);
                histogram(T.sigma(T.autocal==0))
                hold on
                histogram(T.sigma(T.autocal==0),bins);
                title(sprintf('Sigma %.2f +/- %.2f',mean(T.sigma),std(T.sigma)))
                xlabel '\sigma'
                ylabel '#roi'


                nexttile
                % For those where autocal did not work, the tau and amplitude are fixed.
                % Show only the ones that were estimated from the data.
                histogram(T.tau(T.autocal==1))
                hold on
                % Show the default
                line(tpl.parms.deconv.tau*[1 1],ylim,'Color','r','LineWidth',2);
                title(sprintf('Tau %.2f +/- %.2f (def= %.2f)',mean(T.tau(T.autocal==1)),std(T.tau(T.autocal==1)),tpl.parms.deconv.tau))
                xlabel '\tau'
                ylabel '#roi'

                nexttile
                histogram(T.amplitude(T.autocal==1))
                hold on
                line(tpl.parms.deconv.a*[1 1],ylim,'Color','r','LineWidth',2);
                title(sprintf('Amplitude %.2f +/- %.2f (def = %.2f)',mean(T.amplitude(T.autocal==1)),std(T.amplitude(T.autocal==1)),tpl.parms.deconv.a))
                xlabel 'amplitude'
                ylabel '#roi'

                title(layout,sprintf('%s: %d rois from %d sessions in %d subjects',tpl.stag,height(T),numel(unique(T.session_date)),numel(unique(T.subject))));
            end

            if nargout >0
                % Calculate averages per method (stag) , session and
                % subject.
                T = fetchtable(tbl,'*');
                G = groupsummary(T(T.autocal==1,:),["stag" "session_date" "subject"],{@mean,@std},["quality" "amplitude" "sigma" "tau"]);
                nm = strrep(G.Properties.VariableNames,'fun1','mean');
                nm = strrep(nm,'fun2','std');
                G.Properties.VariableNames= nm;
            end
        end

    end

    methods (Access = protected)
        function makeTuples(tbl,key)
            % Select a subset of the roi in the session, run
            % autocalibration and then pick the median of the calibrated 
            % amplitude,sigma, tau and store this in the table for later use 
            % when populating the sbx.Spikes table.
            
            calResults = sbx.mlspike(key,struct([]),calibration =true);
            key.tau = mean([calResults.tau],"omitmissing");
            key.sigma  =mean([calResults.sigma],"omitmissing");
            key.a = mean([calResults.a],"omitmissing");
            key.quality = mean([calResults.quality],"omitmissing");
            key.failed = sum(isnan([calResults.quality]));
            % Store the calibration results             
            insert(tbl, key);

        end
    end
end