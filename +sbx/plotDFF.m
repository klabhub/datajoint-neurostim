function plotDFF(roi,expt,pv)
            % Plot a summary view of dF/F for a set of ROIs in an
            % experiment
            %
            % roi - A sbx.Roi table
            % expt - The ns.Experiment
            % Parm/Value Pairs:
            % baseline - The time window in seconds relative to firstFrame
            %               that defines the  baseline [2 3]
            % respons  - The time window where a response is expected
            %               [0.5 1.5]
            % window    - The time window for which to show dF/F. [ 0 3]
            % maxDFF    - Clamp dF/F values (in %) higher than this [Inf]
            % maxResponse -
            % fetchOptions - Passed to fetch(sbx.Roi), so for instance
            %                   'LIMIT 100' or 'ORDER BY radius'
            % percentile  - Definition of F0 is this percentile of the
            % neuropil corrected fluorescence in the baseline window.
            %                   Set it to 0 to use the mean [8]
            % neuropilFactor - The fraction of hte neuropil to subtract [0.7]
            % shotNoise -Compute shotNoise and show it across the FOV.
            % [false]
            % spikes - Show deconvolved spiking activity instead of
            % Fluoresence.
            % OUTPUT
            %    A figure showing the dF/F for all ROIs, the average dF/F
            %    over time, and a map of the shotnoise
            %
            arguments
                roi (1,1) sbx.Roi {mustHaveRows}
                expt (1,1) ns.Experiment  {mustHaveRows}
                pv.trial   (1,:) double = []
                pv.baseline (1,:) double = [2 3]
                pv.window (1,2) double = [0 3]
                pv.maxdFF (1,1) double = Inf;
                pv.fetchOptions {mustBeText} = ''
                pv.percentile (1,1) double {mustBeInRange(pv.percentile,0,100)} = 8;
                pv.neuropilFactor (1,1) double {mustBeInRange(pv.neuropilFactor,0,1)} = 0.7;
                pv.shotNoise (1,1) logical =false;
                pv.spikes (1,1) logical  =false
            end

            % Extract fluorescence and neuropil
            frame = 1./unique([fetch(sbx.Preprocessed & roi,'framerate').framerate]);
            if pv.spikes
                tF  = get(roi,expt,trial = pv.trial, fetchOptions= pv.fetchOptions, modality='spikes',start = pv.window(1) ,stop=pv.window(2), step = frame,interpolation='nearest');
                tFNeu =[];
                pv.shotNoise =false;
            else
                tF  = get(roi,expt,trial = pv.trial, fetchOptions= pv.fetchOptions, modality='fluorescence',start = pv.window(1) ,stop=pv.window(2), step = frame,interpolation='nearest');
                tFNeu = get(roi,expt,trial = pv.trial,fetchOptions= pv.fetchOptions, modality='neuropil',start = pv.window(1) ,stop=pv.window(2), step = frame,interpolation='nearest');
            end
            % Compute df/f
            if pv.shotNoise
                [dFF,shotNoise] =  sbx.dFOverF(tF,tFNeu,pv.baseline,percentile = pv.percentile,neuropilFactor=pv.neuropilFactor);
            else
                dFF  =  sbx.dFOverF(tF,tFNeu,pv.baseline,percentile = pv.percentile,neuropilFactor=pv.neuropilFactor);
            end
            [~,nrRoi] = size(dFF);


            % Clamp
            dFF(dFF < -abs(pv.maxdFF)) =-abs(pv.maxdFF);
            dFF(dFF > abs(pv.maxdFF)) =abs(pv.maxdFF);

            %% Graphical output
            subplot(2,2,[1 3])
            imagesc(seconds(tF.Time),1:nrRoi,dFF')
            set(gca,'CLim',[0 prctile(dFF(:),95)]);
            xlabel 'Time (s)'
            ylabel (char('ROI',pv.fetchOptions))
            colorbar
            colormap hot
            if pv.shotNoise
                title (sprintf('ShotNoise %.2f +/- %.2f',mean(shotNoise,2,"omitnan"),std(shotNoise,0,2,"omitnan")))
            end
            subplot(2,2,2);
            m = median(dFF,2,"omitnan");
            %e = std(dFF,0,2,"omitnan")./sqrt(sum(~isnan(dFF),2,"omitnan"));
            e = iqr(dFF,2);
            ploterr(tF.Time,m,e)
            hold on
            if ~isempty(pv.baseline)
                patch(seconds([pv.baseline(1) pv.baseline(1) pv.baseline(2) pv.baseline(2)]), [ylim fliplr(ylim)],0.8*ones(1,3),'FaceAlpha',0.1);
            end
            xlabel 'Time (s)'
            ylabel 'dF/F (%)'
            if pv.shotNoise
                subplot(2,2,4)
                if isempty(pv.fetchOptions)
                    roiUsed = roi;
                else
                    roiUsed = roi  & fetch(roi,pv.fetchOptions);
                end
                plotSpatial(roiUsed,color=shotNoise);
            end
            info = fetch(expt,'paradigm');
            sgtitle(sprintf('%s (%s) - %s',info.session_date,info.starttime,info.paradigm))
        end