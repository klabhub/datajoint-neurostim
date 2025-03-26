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
            arguments
                tbl (1,1) sbx.Spikes {mustHaveRows(tbl)}
            end
            sessionFolder = folder(ns.Experiment & tbl);
            planeFolder = fetchtable(sbx.PreprocessedRoi & tbl,'plane');
            filename = fetchtable(tbl,'stag')
        end

    end



    methods (Access=protected)
        function makeTuples(tbl,key)

            assert(exist('fn_structmerge','file'),"The brick repository must be on the path for mlSpike");
            assert(exist('spk_est.m','file'),"The spikes repository must be on the path for mlSpike");
            parms =fetch1(sbx.SpikesParm &key,'parms');
            prep =fetch(sbx.Preprocessed & key,'*');
            %Start parpool if requested
            if parms.nrWorkers>0
                pool = gcp("nocreate");
                if isempty(pool)
                    % Create a pool
                    pool = parpool(parms.nrWorkers);
                else
                    % Use what is available
                    parms.nrWorkers = max(parms.nrWorkers,pool.NumWorkers);
                end
            else
                pool = [];
            end


            sessionPath=unique(folder(ns.Experiment & key));
            dt = 1/prep.framerate;
            parms.deconv.dt = dt;
            roiInPrep = fetch(sbx.PreprocessedRoi & key & parms.restrict,'plane','roi');
            planes = unique([roiInPrep.plane]);
            % Do MLSpike deconvolution
            for p=planes(:)'
                fldr = fullfile(sessionPath,prep.folder,sprintf('plane%d',p));
                % Get the segmented data
                F = ndarrayToArray(py.numpy.load(fullfile(fldr,'F.npy'),allow_pickle=true));
                Fneu = ndarrayToArray(py.numpy.load(fullfile(fldr,'Fneu.npy'),allow_pickle=true));
                roi = [roiInPrep.roi];
                stayRoi = roi; % & plane...
                F =F(stayRoi,:);
                Fneu =Fneu(stayRoi,:);
                signal = F -0.7*Fneu;
                signal = signal';
                [nrSamples,nrRoi] = size(signal);


                %  signal  = signal(1:floor(nrSamples/5),1:2);
                %  [nrSamples,nrRoi] = size(signal);
                %  roi = roi(1:2);
                % %
                spikeCount = nan(nrSamples,nrRoi);
                drift = nan(nrSamples,nrRoi);

                tpl  =struct('roi',num2cell(roi),'nanfrac',nan,'quality',nan,'autocal',nan,'amplitude',nan,'tau',nan,'sigma',nan);
                fprintf('Deconvolving %d channels \n',nrRoi)
                %% Loop for/parfor per channel
                if isempty(pool)
                    for  ch = 1:nrRoi
                        [spikeCount(:,ch),drift(:,ch), quality, nanFrac,cal,autoCal]   = sbx.Spikes.mlSpikeSingleRoi(signal(:,ch),parms);
                        tpl(ch).quality = quality;
                        tpl(ch).nanfrac = nanFrac;
                        tpl(ch).autocal = autoCal;
                        tpl(ch).amplitude = cal.a;
                        tpl(ch).tau = cal.tau;
                        tpl(ch).sigma = cal.finetune.sigma;
                    end
                else
                    parfor  (ch = 1:nrRoi,parms.nrWorkers)
                        dj.conn; % Need to refresh connection in each worker
                        [spikeCount(:,ch),drift(:,ch), quality, nanFrac,cal,autoCal]  =  sbx.Spikes.mlSpikeSingleRoi(signal(:,ch),parms); %#ok<PFOUS>
                        tpl(ch).quality = quality;
                        tpl(ch).nanfrac = nanFrac;
                        tpl(ch).autocal = autoCal;
                        tpl(ch).amplitude = cal.a
                        tpl(ch).tau = cal.tau;
                        tpl(ch).sigma = cal.finetune.sigma;
                    end
                end
                % Save results in a mat file in the plane folder, together
                % with F.npy etc.
                trgFile= fullfile(fldr,[key.stag '.mat']);
                save(trgFile,'spikeCount','drift','-v7.3');
                tpl = mergestruct(key,tpl);
                insert(tbl,tpl);
            end % for plane
        end
    end


    methods (Static)
        %% Function that does the deconvolution for one channel
        % This is called from eitehr a for (nrWorkers=0) or parfor loop
        function [spk,drift, quality, nanFrac,calibratedParms,autoCal] = mlSpikeSingleRoi(signal,parms)
            isNaN = isnan(signal);
            nrSamples = numel(signal);
            signal(isNaN) = 0; % Could remove samples instead or linearly interpolate, but this should be rare (missing F)
            if isfield(parms,'autocalibration') && ~isempty(fieldnames(parms.autocalibration))
                % Do autocalibration
                pax = spk_autocalibration('par'); % Get defaults, then overrule with parms.autocalibrate
                % Copy values from parms.autocalibration struct to pax
                pax = fn_structmerge(pax,parms.autocalibration,'strict','recursive','type');
                pax.mlspikepar = parms.deconv; % Autocalibration calls tps_mlspikes - with these parms
                pax.dt = parms.deconv.dt; % Dont allow overrule by autocalibrate.
                % perform auto-calibration
                [tau,amp,sigma,events] = spk_autocalibration(signal,pax);
                calibratedParms = parms.deconv; % Default from CParm
                calibratedParms.finetune.sigma = sigma;
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
            [thisSpk,fit,drift,parEst] = spk_est(signal,calibratedParms);
            calibratedParms.finetune.sigma= parEst.finetune.sigma; % Estimated (if no autocal)

            %% spk - spiketimes in s
            % Convert back to spike counts at the sample rate of the fluorescence
            bins = (0:nrSamples)*parms.deconv.dt;
            spk = histcounts(thisSpk,bins)';
            %% fit - fit of the F signal at each sample - used to estimate quality

            quality = double(corr(fit(~isNaN),signal(~isNaN),Type="Pearson"));
            nanFrac = mean(isNaN);

        end

    end
end
