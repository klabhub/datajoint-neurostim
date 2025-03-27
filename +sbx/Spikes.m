%{
# Deconvolved Spikes of a ROI in a session, based on a Preprocessed set.
-> sbx.Preprocessed
-> sbx.SpikesParm
roi : smallint
---
quality : float #
nanfrac : float #
autocal : float #
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
%                    'parms',struct('what','mlspikes','prep','gcamp6s'));   % pass these  parms to the fun
%
% insertIfNew(ns.CParm, spksPrep);
% and then call
% populate(ns.C,'ctag="spikes"')
%
%

classdef Spikes < dj.Computed
    properties  (Dependent)
        file
    end

    methods
        function v = get.file(tbl)
            % Returns the name of the file where the MLSpike results are
            % saved. Note that this looks only at the first row in tbl; all
            % subsequent rows are ignored
            arguments
                tbl (1,1) sbx.Spikes {mustHaveRows(tbl)}
            end
            tpl =fetch(tbl,'LIMIT 1');
            fldr = unique(folder(ns.Session & tpl));
            plane= fetchn(sbx.PreprocessedRoi & tpl,'plane');
            planeFldr = sprintf('plane%01d',unique(plane));
            v = fullfile(fldr,planeFldr,[tpl.stag '.mat']);
        end
    end
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
                        
            warning('off','backtrace');
            parms =fetch1(sbx.SpikesParm &key,'parms');  % Parameters for mlspike
            prep =fetch(sbx.Preprocessed & key,'*');     % Preprocessed data set
            
            pool =nsParPool();

            sessionPath=unique(folder(ns.Experiment & key));
            dt = 1/prep.framerate;
            parms.deconv.dt = dt;
            % Restrict rois as specified in parms.restrict
            roiInPrep = fetch(sbx.PreprocessedRoi & key & parms.restrict,'plane','roi');
            planes = unique([roiInPrep.plane]);
            % Do MLSpike deconvolution per roi
            for p=planes(:)'
                fldr = fullfile(sessionPath,prep.folder,sprintf('plane%d',p));
                trgFile = fullfile(fldr,[key.stag '.mat']);
                if exist(trgFile,"file")
                    fprintf("%s already exists, reading from file.\n",trgFile);
                    secs = load(trgFile,'info');
                    info = secs.info;
                else
                    % Get the segmented data
                    fprintf('Loading segmentation from %s.\n',fldr)
                    F = ndarrayToArray(py.numpy.load(fullfile(fldr,'F.npy'),allow_pickle=true));
                    Fneu = ndarrayToArray(py.numpy.load(fullfile(fldr,'Fneu.npy'),allow_pickle=true));
                    roi = [roiInPrep.roi];                    
                    F =F(roi,:);
                    Fneu =Fneu(roi,:);
                    signal = F -0.7*Fneu; % Subtract neuropil
                    signal = signal'; % Rois as columns for parfor
                    [nrSamples,nrRoi] = size(signal);
                    spikeCount = nan(nrSamples,nrRoi);
                    drift = nan(nrSamples,nrRoi);
                    info  =struct('roi',num2cell(roi),'nanfrac',nan,'quality',nan,'autocal',nan,'amplitude',nan,'tau',nan,'sigma',nan);
                    fprintf('Deconvolving %d channels \n',nrRoi)
                    %% Loop for/parfor per channel
                    dq = parallel.pool.DataQueue;
                    counter = 0;     totalDuration =minutes(0);                                        
                    afterEach(dq, @(x) updateMessage(x));
                    tStart = tic;
                    if isempty(pool)                        
                        for  ch = 1:nrRoi
                            tic;
                            send(dq,{ch,false,0});
                            [spikeCount(:,ch),drift(:,ch), quality, nanFrac,cal,autoCal]   = sbx.Spikes.mlSpikeSingleRoi(signal(:,ch),parms);
                            % Store autocalibration results and other info
                            info(ch).quality = quality;
                            info(ch).nanfrac = nanFrac;
                            info(ch).autocal = autoCal;
                            info(ch).amplitude = cal.a;
                            info(ch).tau = cal.tau;
                            info(ch).sigma = cal.finetune.sigma;
                            send(dq,{ch,true,seconds(toc)});
                        end
                    else 
                        parfor  (ch = 1:nrRoi)
                            dj.conn; % Need to refresh connection in each worker
                            warning('off','backtrace'); % Needs to be set on each worker
                            tic
                            send(dq,{ch,false,0})
                            [spikeCount(:,ch),drift(:,ch), quality, nanFrac,cal,autoCal]  =  sbx.Spikes.mlSpikeSingleRoi(signal(:,ch),parms); %#ok<PFOUS>
                            info(ch).quality = quality;
                            info(ch).nanfrac = nanFrac;
                            info(ch).autocal = autoCal;
                            info(ch).amplitude = cal.a
                            info(ch).tau = cal.tau;
                            info(ch).sigma = cal.finetune.sigma;  
                            send(dq,{ch,true, seconds(toc)})
                        end
                    end
                    % Save results in a mat file in the plane folder, together
                    % with F.npy etc. This file will be read by sbx.read to add
                    % the spike as continuous data in ns.C
                    info = mergestruct(key,info);
                    save(trgFile,'spikeCount','drift','info','-v7.3');
                end
                % Store meta data in the table.
                insert(tbl,info);
            end % for plane      
            warning('on','backtrace');
            function updateMessage(x)
                [channel,done,thisDuration] =deal(x{:});
                if done
                    counter= counter+1;                                               
                    fprintf("Deconvolution complete (%d out of %d : %s, cumulative %s) \n",counter,nrRoi,thisDuration,minutes(toc(tStart)));
                else
                    fprintf("Starting channel #%d\n",channel);
                end
            end
        end
    end

    methods (Static)
        %% Function that does the deconvolution for one channel
        % This is called from either a for or parfor loop; set NS_PARFOR to
        % the number of workers or create a parpool manually to use parfor.
        function [spk,drift, quality, nanFrac,calibratedParms,autoCal] = mlSpikeSingleRoi(signal,parms)
            isNaN = isnan(signal);
            nrSamples = numel(signal);
            signal(isNaN) = 0; % Could remove samples instead or linearly interpolate, but this should be rare (missing F)
            if isfield(parms,'autocalibration') && ~isempty(fieldnames(parms.autocalibration))
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
            else
                % No autocalibration : use parms as specified in CParm
                calibratedParms = parms.deconv;
                autoCal = false;
            end
            % Do the deconvolution
            [thisSpk,fit,drift,parEst] = spk_est(signal,calibratedParms);
            calibratedParms.finetune.sigma= parEst.finetune.sigma; % Estimated (if no autocal)

            % Convert back to spike counts at the sample rate of the fluorescence
            bins = (0:nrSamples)*parms.deconv.dt;
            spk = histcounts(thisSpk,bins)';
            % fit - fit of the F signal at each sample - used to estimate quality
            quality = double(corr(fit(~isNaN),signal(~isNaN),Type="Pearson"));
            nanFrac = mean(isNaN);
        end

    end
end
