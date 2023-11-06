function plot(roi,expt, condition,pv)
            % Plot time courses, spectra. or tuning.
            %
            % Each roi in the table will be shown as a
            % separate tile, each condition a line in the plot.  Time
            % courses are scaled to the 99th percentile across all
            % responses and de-meaned per condition. Hence, the mean
            % response is lost, but the relative response modulation in
            % each condition is maintained.
            %
            % roi - sbx.Roi table
            % expt - A single ns.Experiment (tuple or table)
            % condition - Specify how trials should be grouped into conditions:
            %               []  - Pool over all trials
            %               A ns.Condition table - pool per condition
            %               A vector of trials - Pool over only these trials
            %               A cell array with vectors of trials. Pool over
            %               each set of trials
            %
            % 'fun' - By default this function visualizes the mean across trials
            %       together with shading reflecting the standard error. To use something else,
            %       pass a function that, when passed a matrix with [nrTimePoints nrTrials] , returns one
            %       value and an error bar for each row.
            %
            % 'name'  Name of the conditions
            % 'start' - Start time in seconds
            % 'step'  - Step time in seconds
            % 'stop' - Stop time in seconds
            % 'interpolation' - Interpolation method ['linear']
            % 'crossTrial ' - Allow start/stop to cross to the
            %               previous/next trial.
            % 'mode'  ["TIMECOURSE"], TUNING, SPECTRUM ,"RASTER"
            % 'perTrial'  -Show individual trials [false]
            % 'prctileMax'  Percentile that is used to scale responses for
            %               visualization [95]
            % Spectrum Options
            % 'evoked' Set  to true to show evoked power instead of total
            % power.
            % Tuning options
            %   pv.x    - The independent variable for each condition
            %  pv.polar . Set to true to indicate that pv.x is in degrees.

            arguments
                roi (1,1) sbx.Roi {mustHaveRows}
                expt (1,1) ns.Experiment {mustHaveRows}
                condition = []
                pv.fun (1,1) = @(x)(deal(mean(x,2,"omitnan"),std(x,0,2,"omitnan")./sqrt(sum(~isnan(x),2))));
                pv.name {mustBeText} = {}
                pv.start (1,1) double = 0
                pv.stop (1,:) double =  3
                pv.step  (1,1) double = 1/15.5;
                pv.interpolation {mustBeText} = 'linear';
                pv.modality {mustBeText} = 'spikes';
                pv.averageRoi (1,1)  logical = false;
                pv.mode (1,1) {mustBeTextScalar,mustBeMember(pv.mode,["COHERENCE", "RASTER", "TIMECOURSE","EVOKED","TOTAL", "TUNING"])} = "TIMECOURSE"
                pv.crossTrial (1,1) logical = false;
                pv.fetchOptions {mustBeText} = ''
                pv.perTrial (1,1) logical = false;
                pv.prctileMax (1,1) double {mustBeInRange(pv.prctileMax,0,100)} = 95;
                % Layout
                pv.compact = false;
                % Spectrum options
                pv.evoked (1,1) logical = false;
                pv.options cell = {}; % Cell array of parameter value pairs passed to pspectrum
                % Tuning options
                pv.x (1,:) {mustBeNumeric} =[]
                pv.polar (1,1) logical = false;
            end

            [trPerCondition,names] = trialsPerCondition(condition);
            if isempty(pv.name)
                pv.name = names;
            end

            nrConditions = numel(trPerCondition);
            if isempty(pv.name)
                pv.name= "Condition " + string(1:nrConditions);
            end

            if isempty(pv.fetchOptions)
                roiTpls = fetch(roi);
            else
                roiTpls = fetch(roi,pv.fetchOptions);
            end
            nrRois = numel(roiTpls);
            if pv.averageRoi || pv.mode=="COHERENCE"
                nrRois =1;
            end
            layout = tiledlayout('flow');
            if pv.compact
                layout.Padding ="tight";
            end
            % Loop over rois
            for roiCntr = 1:nrRois
                m = [];
                e = [];
                allTime = [];
                perTrial =cell(1,nrConditions);
                for c= 1:nrConditions
                    stop = pv.stop(min(c,numel(pv.stop)));
                    % Loop over conditions
                    if c==1
                        nexttile;
                    end
                    if pv.averageRoi || pv.mode=="COHERENCE"
                        % Get all rois
                        [time,y] = get(roi ,expt,fetchOptions = pv.fetchOptions,crossTrial =pv.crossTrial, trial=trPerCondition{c},modality = pv.modality,start=pv.start,stop=stop,step=pv.step,interpolation =pv.interpolation);
                        if pv.averageRoi
                            % Average
                            y = mean(y,2,"omitnan"); % Average over rois
                        end
                    else
                        %% Loop over roi, one tile per roi
                        [time,y] = get(roi &roiTpls(roiCntr) ,expt,crossTrial =pv.crossTrial,trial=trPerCondition{c},modality = pv.modality,start=pv.start,stop=stop,step=pv.step,interpolation =pv.interpolation);
                    end

                    if isempty(y);continue;end

                    perTrial{c} = y;
                    switch upper(pv.mode)
                        case {"TOTAL","EVOKED"}
                            y = y-mean(y,1,"omitnan");
                            if upper(pv.mode) =="EVOKED"
                                y = mean(y,2,"omitnan");
                            end
                            y(isnan(y)) =0;
                            [pwr,freq] = pspectrum(y,time,'power',pv.options{:});
                            [thisM,thisE] = pv.fun(pwr); % Average over trials
                        case "TUNING"
                            window = isbetween(seconds(time),seconds(pv.start),seconds(stop));
                            meanResponseInWindow = squeeze(mean(y(window,:),1,"omitnan")); % Average over window
                            [thisM,thisE] = pv.fun(meanResponseInWindow);
                        case {"TIMECOURSE", "RASTER"}
                            % Average over trials  in the condition
                            [thisM,thisE] = pv.fun(y);
                        case "COHERENCE"
                            y = y- mean(y,1,"omitnan"); % Remove mean
                            y(isnan(y)) = 0; % Remove nans
                            for tr =  1:size(y,2)
                                [thisC(:,:,:,tr),phi,S12,freq] = cohmatrixc(squeeze(y(:,tr,:)),struct('tapers',[3 5],'pad',0,'Fs',1./pv.step));
                            end
                            thisM = mean(thisC,4);
                    end
                    m = catpad(m,thisM);  % Cat as next column allow different rows (padded with NaN at the end)
                    e  =catpad(e,thisE);
                    if numel(time)>numel(allTime)
                        allTime  =time;
                    end
                end

                axes(layout.Children(1)); %#ok<LAXES>

                %% Visualize
                switch upper(pv.mode)
                    case {"TOTAL","EVOKED"}
                        h = ploterr(freq,squeeze(m),squeeze(e),'linewidth',2,'ShadingAlpha',0.5);
                        xlabel 'Frequency (Hz)'
                        ylabel (pv.mode + ' Power')
                        h =legend(h,pv.name);
                        h.Interpreter  = 'None';

                    case "TUNING"
                        if isempty(pv.x)
                            pv.x = (1:nrConditions)';
                        end
                        [pv.x,ix] =sort(pv.x);
                        m = m(ix)';
                        e = e(ix)';
                        pv.name = pv.name(ix);
                        if pv.polar
                            %Assume uX is degrees. Show in polar coordinates,
                            %connect lines around the circle.
                            x= deg2rad(pv.x)';
                            x= [x;x(1)];
                            y = [m;m(1)];
                            e = [e;e(1)];
                            polarplot(x,y);
                            hold on
                            polarplot(x,y-e,'k:');
                            polarplot(x,y+e,'k:');
                        else
                            ploterr(pv.x,m,e)
                            set(gca,'xTick',pv.x,'xTickLabel',pv.name)
                        end

                    case "RASTER"
                        grandMax = max(cellfun(@(x) prctile(x(:),pv.prctileMax ),perTrial));
                        grandMin = min(cellfun(@(x) prctile(x(:),100-pv.prctileMax ),perTrial));
                        nrTime = numel(allTime);
                        cmap = [hot(255);0 0 1];
                        I =[];
                        for c=1:nrConditions
                            I = cat(2,I,perTrial{c},nan(nrTime,1));
                        end



                        I = ((I-grandMin)./(grandMax-grandMin))';
                        % Clamp
                        I(I<0) = 0;
                        I(I>1) = 1;
                        I= round(I*255);
                        I(isnan(I))=256;
                        nrTrials = size(I,1);
                        image(allTime,1:nrTrials, I,'CDataMapping','direct');
                        colormap(cmap)
                        title(['ROI #' num2str(roiTpls(roiCntr).roi)]);
                        nrTrialsPerCondition = cellfun(@(x) size(x,2),perTrial);
                        leftEdge = [0 cumsum(nrTrialsPerCondition(1:end-1))];
                        middleOfCondition = leftEdge+nrTrialsPerCondition./2;
                        set(gca,'yTick',middleOfCondition,'yTickLabel',pv.name)

                        if pv.compact
                            set(gca,'XTick',[]);
                        else
                            xlabel 'Time (s)'
                            ylabel('Conditions')
                            h = colorbar;
                            set(h,'YTick',0:50:250,'YTickLabel',round(grandMin +(0:50:250)*(grandMax-grandMin)/255) )
                            ylabel(h,'Response')
                        end

                    case "TIMECOURSE"



                        %% TimeCourse
                        nrTime = numel(allTime);
                        % Scale each condition to the grandMax
                        if pv.perTrial
                            grandMax = max(cellfun(@(x) prctile(x(:),pv.prctileMax ),perTrial));
                            grandMin = min(cellfun(@(x) prctile(x(:),100-pv.prctileMax ),perTrial));
                        else
                            grandMax = prctile(abs(m(:)),pv.prctileMax );
                            grandMin = prctile(abs(m(:)),100-pv.prctileMax );
                        end
                        m = (m-grandMin)./(grandMax-grandMin);
                        e = e./(grandMax-grandMin);
                        % Add the conditionNr so that each m column has a mean of
                        % conditionNr and can be plotted on the same axis, with
                        % conditions discplaced vertically from each other.
                        m = m + repmat(1:nrConditions,[nrTime 1]);
                        [h,hErr] = ploterr(allTime,m,e,'linewidth',2,'ShadingAlpha',0.5);
                        hold on
                        % Show "zero" line
                        hh = plot(allTime,repmat(1:nrConditions,[nrTime 1]),'LineWidth',0.5);
                        [hh.Color] =deal(h.Color);
                        ylim([1 nrConditions+1])
                        set(gca,'yTick',1:nrConditions,'yTickLabel',pv.name)
                        xlabel 'Time (s)'
                        ylabel 'Response per condition'
                        if pv.perTrial
                            [hh.Color] = deal([0 0 0]);
                            [h.Color] = deal([0 0 0]);
                            [hErr.FaceColor] = deal([0 0 0]);
                            colorOrder = get(gca,'ColorOrder');
                            for c=1:nrConditions
                                plot(allTime,c+ perTrial{c}./(grandMax-grandMin),'Color',colorOrder(mod(c-1,size(colorOrder,1))+1,:),'LineWidth',.5)
                            end

                        end
                end
            end
        end