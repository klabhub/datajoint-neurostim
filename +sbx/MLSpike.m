%{
#  Autocalibrated parameter settings for MLSpike, based on a Preprocessed set and initial values in MLSpikeParm.
-> sbx.Preprocessed   # Preprocessed set used to autocalibrate 
-> sbx.SpikesParm     # Non-autocalibrated parameters for the MLSPike algorithm
---
quality : float #   Average quality
nanfrac : float #   
amplitude : float
sigma : float
tau : float
%}
%
% Perform spike deconvolution on a Calcium fluorescence signal using the MLspike algorithm
% Deneux, T., Kaszas, A., Szalay, G., Katona, G., Lakner, T., Grinvald, A., Rozsa, B., & Vanzetta, I. (2016). Accurate spike estimation from noisy calcium signals for ultrafast three-dimensional imaging of large neuronal populations in vivo. Nature Communications 2016 7:1, 7(1), 1–17. https://doi.org/10.1038/ncomms12190
% Relies on the spikes and brick repositories: https://github.com/MLspike
%
%  See sbx.SpikesParm for the definition of the parms struct
% After creating a sbx.Preprocessed set and defining the SpikesParm entry,
% populate this table. The SpikesParm has a restrict field to allow limiting deconvolution
%  to cells (e.g. 'pcell>0.75')
%
% Populating this table will perform the deconvolution and store the
% results in a file named after the stag, in the same folder as the
% (suite2p) segmentation results (.file property of this class)
% The table only contains meta data, not the actual spikes. To use the spikes,
% populate a ns.C table with
% spksPrep = struct('ctag','spikes',...       % Name of this C
%                   'description','ML Spike deconvolved spikes ',...
%                    'extension','.sbx',...  % Files to process
%                    'fun','sbx.read',...      % Use this function to read the spikes files
%                    'parms',struct('what','mlspikes','prep','gcamp6s'));  
% where 'mlspikes'  matches the stag used to perform the deconvolution
% (i.e. the row in the sbx.Spikes table)
%
% insertIfNew(ns.CParm, spksPrep);
% and then call
% populate(ns.C,'ctag="spikes"')
%
%

classdef  MLSpike< dj.Computed
 
    methods
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

    methods (Access=protected)
        function makeTuples(tbl,key)
            % Confirm that the mlspike toolbox is on the path
            assert(exist('fn_structmerge','file'),"The brick repository must be on the path for mlSpike");
            assert(exist('spk_est.m','file'),"The spikes repository must be on the path for mlSpike");
            setupPython;
            warning('off','backtrace');
            parms =fetch1(sbx.MLSpikeParm &key,'parms');  % Parameters for mlspike
            prep =fetch(sbx.Preprocessed & key,'*');     % Preprocessed data set

            candidateRoi = sbx.PreprocessedRoi & key & parms.restrict;
            nrRoi = count(candidateRoi);
            nrRoiToUse = round(max(0.05*nrRoi,min(nrRoi,50)));
            channels = fetchn(candidateRoi,"ROI");
            channels = channels(randperm(nrRoi));
            channels = channels(1:nrRoiToUse);


            dt = 1/prep.framerate;
            parms.deconv.dt = dt;
            tStart = tic;
            for  ch = channels
                tic;
                F  = fetchn(ns.C & 'ctag="fluorescence"' & key & struct('channel',ch),'signal');


                      isNaN = isnan(F);
            nrSamples = numel(F);
            F(isNaN) = 0; % Could remove samples instead or linearly interpolate, but this should be rare (missing F)
            % Do autocalibration
            pax = spk_autocalibration('par'); % Get defaults, then overrule with parms.autocalibrate
            % Copy values from parms.autocalibration struct to pax
            pax = fn_structmerge(pax,parms.autocalibration,'strict','recursive','type');
                % For consistency of autocalibration and deconvolution;
            % copy deconv to mlspike par, the autocalibration function
            % calls tps_mlspikes with these parms
            pax.mlspikepar = parms.deconv;
            pax.dt = parms.deconv.dt; % Dont allow overrule by autocalibrate.
            % perform auto-calibration
            [tau,amp,sigma,events] = spk_autocalibration(signal,pax);
            calibratedParms = parms.deconv; % Default from CParm
            calibratedParms.finetune.sigma = sigma; % Always estimated
                % If events is empty, calibration was not possible, use defaults.
                if isempty(events)
                    autoCal = false; % Failed
                else
                    calibratedParms.tau = tau;
                    calibratedParms.a = amp;
                    autoCal = true;
                end
            
                % Do the deconvolution
            [thisSpk,fit,drift,parEst] = spk_est(signal,calibratedParms);
            calibratedParms.finetune.sigma= parEst.finetune.sigma; % Estimated (if no autocal)

            % Convert back to spike counts at the sample rate of the fluorescence
            bins = (0:nrSamples)*parms.deconv.dt;
            spikeCount = histcounts(thisSpk,bins)';
            % fit - fit of the F signal at each sample - used to estimate quality
            quality = double(corr(fit(~isNaN),signal(~isNaN),Type="Pearson"));
            nanFrac = mean(isNaN);

            % Store autocalibration results and other info
            info.quality = quality;
            info.nanfrac = nanFrac;
            info.autocal = autoCal;
            info.amplitude = calibratedParms.a;
            info.tau = calibratedParms.tau;
            info.sigma = calibratedParms.finetune.sigma;
            info.roi = roi;

      
                sbx.MLSpike.mlSpikeSingleRoi(F,parms,roiToProcess(ch),tempFolder);
                            

                    % Combine individual files into a single trgFile
                    fprintf('Combining individual ROI files into %s.\n', trgFile);
                    spikeCount = nan(nrSamples,nrRoi);
                    drift = nan(nrSamples,nrRoi);
                    info  =struct('roi',num2cell(roi),'nanfrac',nan,'quality',nan,'autocal',nan,'amplitude',nan,'tau',nan,'sigma',nan);
                    for ch  = 1:numel(roi)
                        tempData = load(fullfile(tempFolder,sprintf('roi_%d.mat', roi(ch))));
                        spikeCount(:,ch) =tempData.spikeCount;
                        drift(:,ch) = tempData.drift;
                        info(ch) =  tempData.info;
                    end

                    % Save combined data to trgFile
                    save(trgFile, 'spikeCount', 'drift', 'info', '-v7.3');
                    fprintf('Saved combined data to %s.\n', trgFile);
                end
                
                % Store the results in the DJ database
                info = mergestruct(key,info);
                insert(tbl,info);
end % for plane


            warning('on','backtrace');
            
        end
    end

    
end
