%{
# ROI in a session with is Fluorescence and spiking activity.
-> sbx.Preprocessed
roi : smallint    #  id within this segmentation and plane
plane : smallint #  plane
---
pcell =0 : float # Probability that the ROI is a neuron
fluorescence   : longblob     # Fluorescence trace
neuropil: longblob   # Neuropil Fluorescence trace (scatter)            
spikes: longblob   # Deconvolved spiking activity
meanrate : float  # Mean number of spikes/second across the session
stdrate  :  float  # Stdev of the spikes/second across the session
x        : Decimal(4,0) # x-Pixel location 
y        : Decimal(4,0) # y-Pixel location 
radius   : Decimal(4,1) # Radius of cell in micron
aspect   : Decimal(4,1)  # Aspect raio of major/minor axis of a 2D Gaussian fit.
compact  : Decimal(4,2)  # How compact the ROI is ( 1 is a disk, >1 means less compact)
%}
classdef Roi < dj.Imported

    methods (Access = public)
        function [hScatter,ax,axImg] = plot(roi,pv)
            % Show properties of ROIs in a spatial layout matching that
            % used in suite2p.
            %
            % INPUT
            % roi  = (subset of) sbx.Roi
            % Parameter/Value pairs:
            % sz = Size to use for each ROI. By default, the estimated physical size
            %       of the roi is used. But by passing some other property (e.g. tuning per ROI),
            %       this can be visualized. Maximum size is 500 points (surface area).
            % szLabel  - Label to use in the datatip.
            % color = Color to use for each ROI. By default, the z-scored
            %           spiking activity across the session is used.
            % colorLabel = Label to use in the datatip.
            % clim - Color limits to use
            % showImg  - Set to true to show the mean Ca image as a
            %               background [true]
            % colormap - Colormap to use  [hot]
            %
            % OUTPUT
            %    hScatter - Handle to the scatter object
            %   ax - Axes showing the cells
            %  axImg - Axes showing the backgroun image.
            % EXAMPLE
            %
            %
            % plot(sbx.Roi); % Encode all cell radius (size) and z-scored session activity (color).
            % z = computed properties of all rois
            % plot(sbx.Roi,color =z,clim = [0 5]);  % Encode z by the color, crop colors at 5
            arguments
                roi (1,1) sbx.Roi
                pv.sz      = []
                pv.szLabel = ''
                pv.color   = []
                pv.colorLabel = ''
                pv.clim     =[];
                pv.showImg = true;
                pv.colormap = hot;
            end
            MAXPOINTS = 500; % S
            [x,y,radius,m,sd] = fetchn(roi,'x','y','radius','meanrate','stdrate');
            micPerPix = sbx.micronPerPixel(roi);  % Scaling
            mixPerPixR = sqrt(sum(micPerPix.^2));
            % Setup size
            if isempty(pv.sz)
                % Use the physical size as the size of the cells
                pv.sz =min(pi*(radius*mixPerPixR).^2,MAXPOINTS);
                pv.szLabel = 'size (\mum^2)';
            end
            zeroSize = pv.sz==0;
            if any(zeroSize)
                pv.sz(zeroSize) = eps;
            end
            % Setup color
            if isempty(pv.color)
                % Use the z-scored rate as the color of the cells
                z = m./sd;
                pv.color = z;
                pv.colorLabel = 'rate (Z)';
            end
            colormap(pv.colormap)
            if pv.showImg
                % Show mean image as background
                meanImg = fetch1(sbx.Preprocessed & roi,'img','LIMIT 1');
                cl = [min(meanImg,[],"all") max(meanImg,[],"all")];
                RI = imref2d(size(meanImg),micPerPix(1),micPerPix(2));
                imshow(meanImg,RI,'DisplayRange',cl,'colormap',gray);
                axImg =gca;
                % Another set of axes for the scatter (uses a different
                % colormap)
                ax = axes('position',get(axImg,'position'),'Color','none');
            else
                ax =gca;
            end

            % The data
            hScatter = scatter(ax,y*micPerPix(2),x*micPerPix(1),pv.sz(:), pv.color(:),'filled');
            xlabel 'Y (\mum)'
            ylabel 'X (\mum)'
            % Match Suite2p layout
            set(ax,'XDir','normal','YDir','reverse');
            if ~isempty(pv.clim)
                clim(ax,pv.clim)
            end
            ax.Color = "none";
            h = colorbar;
            ylabel(h,pv.colorLabel);
            if exist('axImg','var')
                set(axImg,'position',get(ax,'position'));
                ax.XLim = axImg.XLim;
                ax.YLim = axImg.YLim;
                ax.Visible ='off';
            end

            % Data tips
            if ~isempty(pv.szLabel)
                hScatter.DataTipTemplate.DataTipRows(3).Label = pv.szLabel;
            end
            if ~isempty(pv.colorLabel)
                hScatter.DataTipTemplate.DataTipRows(4).Label = pv.colorLabel;
            end

        end

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
            % maxDFF    - Clamp dF/F values (in %) higher than this [100]
            % fetchOptions - Passed to fetch(sbx.Roi), so for instance
            %                   'LIMIT 100' or 'ORDER BY radius'
            % percentile  - Definition of F is this percentile of the
            % neuropil corrected fluorescence in the baseline window. [8]
            % OUTPUT
            %    A figure showing the dF/F for all ROIs, the average dF/F
            %    over time, and a map of the shotnoise
            %
            arguments
                roi (1,1) sbx.Roi
                expt (1,1) ns.Experiment
                pv.trial   =[]
                pv.baseline = [2 3]
                pv.window = [0 3]
                pv.response = [0.5 1.5]
                pv.maxdFF= Inf;
                pv.maxResponse = Inf
                pv.fetchOptions = ''
                pv.percentile = 8;
                pv.neuropilFactor = 0.7;
                pv.shotNoise =false;
                pv.spikes =false
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
            % Remove outliers based on expected response
            if ~isinf(pv.maxResponse)
                inResponse = tF.Time >= seconds(pv.response(1)) & tF.Time < seconds(pv.response(2));
                mResponse = mean(dFF(inResponse,:),1,"omitnan");
                out = mResponse > pv.maxResponse;
                dFF(:,out) =[];
                nrRoi = nrRoi-sum(out);
                fprintf('Removed %d rois based on their response\n',sum(out));
            end

            % Clamp
            dFF(dFF < -abs(pv.maxdFF)) =-abs(pv.maxdFF);
            dFF(dFF > abs(pv.maxdFF)) =abs(pv.maxdFF);

            %% Graphical output
            subplot(2,2,[1 3])
            imagesc(seconds(tF.Time),1:nrRoi,dFF')
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
            patch(seconds([pv.baseline(1) pv.baseline(1) pv.baseline(2) pv.baseline(2)]), [ylim fliplr(ylim)],0.8*ones(1,3),'FaceAlpha',0.1);
            xlabel 'Time (s)'
            ylabel 'dF/F (%)'
            if pv.shotNoise
                subplot(2,2,4)
                if isempty(pv.fetchOptions)
                    roiUsed = roi;
                else
                    roiUsed = roi  & fetch(roi,pv.fetchOptions);
                end
                plot(roiUsed,color=shotNoise);
            end
            info = fetch(expt,'paradigm');
            sgtitle(sprintf('%s (%s) - %s',info.session_date,info.starttime,info.paradigm))
        end

        function plotTimeCourse(roi,expt, condition,pv)
            % Plot time courses. Each roi in the table will be shown as a
            % separate tile, each condition a line in the plot.  Time
            % courses are scaled to the 99th percentile across all
            % responses and de-meaned per condition. Hence, the mean
            % response is lost, but the relative response modulation in
            % each condition is maintained.
            %
            % roi - sbx.Roi table
            % expt - A single ns.Experiment (tuple or table)
            % conditions - Specify how trials should be grouped into conditions:
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
            %
            % 'mode'  ["TIMECOURSE"], TUNING, SPECTRUM
            % Spectrum Options
            % 'evoked' Set  to true to show evoked power instead of total
            % power.
            % Tuning options
            %   pv.x    - The independent variable for each condition
            %  pv.polar . Set to true to indicate that pv.x is in degrees.

            arguments
                roi (1,1) sbx.Roi
                expt (1,1) ns.Experiment
                condition = []
                pv.fun (1,1) = @(x)(deal(mean(x,2,"omitnan"),std(x,0,2,"omitnan")./sqrt(sum(~isnan(x),2))));
                pv.name= {}
                pv.start = 0
                pv.stop =  3
                pv.step = 0.1
                pv.interpolation = 'linear';
                pv.modality = 'spikes';

                pv.mode (1,1) {mustBeTextScalar,mustBeMember(pv.mode,["TIMECOURSE","SPECTRUM", "TUNING"])} = "TIMECOURSE"

                % Spectrum options
                pv.evoked (1,1) logical = false;
                % Tuning options
                pv.x (1,:) {mustBeNumeric} =[]
                pv.polar (1,1) logical = false;
            end

            %% Determine trial->condition mapping
            if isempty(condition)
                trialsPerCondition = {[]};% All trials
            elseif isa(condition,'ns.Condition')
                % Collect trials from ns.Condition table
                trialsPerCondition = cell(1,count(condition));
                cCntr = 0;
                for c=fetch(condition)'
                    cCntr = cCntr+1;
                    trialsPerCondition{cCntr} = [fetch(ns.ConditionTrial & c,'trial').trial];
                    pv.name{cCntr} = c.name;
                end
            elseif isnumeric(condition)
                trialsPerCondition = {condition};
            elseif iscell(condition)
                trialsPerCondition = condition;
            else
                error('Condition must be empty, an ns.Condition table, a list of trials, or a cell array with lists of trials. Not a %s',class(condition))
            end
            nrConditions = numel(trialsPerCondition);
            if isempty(pv.name)
                pv.name= "Condition " + string(1:nrConditions);
            end

            %% Loop over roi, one tile per roi
            for thisRoi =fetch(roi)'
                nexttile;
                for c= 1:nrConditions
                    % Loop over conditions
                    [time,y] = get(roi & thisRoi,expt,trial=trialsPerCondition{c},modality = pv.modality,start=pv.start,stop=pv.stop,step=pv.step,interpolation =pv.interpolation);
                    if isempty(y);continue;end
                    y = reshape(squeeze(y),size(y,1),[]);
                    switch upper(pv.mode)
                        case "SPECTRUM"
                            y = y-mean(y,1,"omitnan");
                            if pv.evoked
                                y = mean(y,2,"omitnan");
                            end
                            y(isnan(y)) =0;
                            [pwr(:,:,c),freq] = pspectrum(y,time,'power'); %#ok<AGROW>
                        case "TUNING"
                            window = isbetween(seconds(time),seconds(pv.start),seconds(pv.stop));
                            meanResponseInWindow = squeeze(mean(y(window,:),1,"omitnan")); % Average over window
                            [m(c,:),e(c,:)] = pv.fun(meanResponseInWindow); %#ok<AGROW>

                        case "TIMECOURSE"
                            % Average over trials  in the condition
                            [m(:,c),e(:,c)] = pv.fun(y); %#ok<AGROW>
                    end
                end

                %% Visualize
                switch upper(pv.mode)
                    case "SPECTRUM"
                        [m,e] = pv.fun(pwr);
                        maxM = max(m);
                        m = m./maxM;
                        e = e./maxM;
                        h = ploterr(freq,squeeze(m),squeeze(e),'linewidth',2,'ShadingAlpha',0.5);
                        xlabel 'Frequency (Hz)'
                        ylabel 'Power'
                        h =legend(h,pv.name);
                        h.Interpreter  = 'None';
                    case "TUNING"
                        if isempty(pv.x)
                            pv.x = (1:nrConditions)';
                        end
                        [pv.x,ix] =sort(pv.x);
                        m = m(ix,:);
                        e = e(ix,:);
                        pv.name = pv.name(ix);
                        if pv.polar
                            %Assume uX is degrees. Show in polar coordinates,
                            %connect lines around the circle.
                            x= deg2rad(pv.x)';
                            x= [x;x(1)];
                            y = [m;m(1,:)];
                            e = [e;e(1,:)];
                            polarplot(x,y);
                            hold on
                            polarplot(x,y-e,'k:');
                            polarplot(x,y+e,'k:');
                        else
                            ploterr(pv.x,m,e)
                            set(gca,'xTick',pv.x,'xTickLabel',pv.name)
                        end

                    case "TIMECOURSE"
                        %% TimeCourse

                        nrTime = numel(time);
                        % Remove mean per condition
                        m = m - mean(m,1,"omitnan");
                        % Scale each condition to the grandMax
                        grandMax = prctile(abs(m(:)),99);
                        m = m./grandMax;
                        e = e./grandMax;
                        % Add the conditionNr so that each m column has a mean of
                        % conditionNr and can be plotted on the same axis, with
                        % conditions discplaced vertically from each other.
                        m = m + repmat(1:nrConditions,[nrTime 1]);
                        h = ploterr(time,m,e,'linewidth',2,'ShadingAlpha',0.5);
                        hold on
                        % Show mean line
                        hh = plot(time,repmat(1:nrConditions,[nrTime 1]),'LineWidth',0.5);
                        [hh.Color] =deal(h.Color);
                        ylim([0 nrConditions+1])
                        set(gca,'yTick',1:nrConditions,'yTickLabel',pv.name)
                        
                end
            end
        end


        function [varargout] = get(roi,expt,pv)
            % Function to retrieve trial-start aligned activity data per
            % ROI and Experiment.
            % roi  - Table o sbx.Roi to use
            % expt - Table of ns.Experiment to use
            % Optional Parameter/Value pairs
            % modality - 'spikes','fluorescence','neuropil'
            % trial   - Which trials to extract
            % start  - First time point to extract (relatve to first frame
            %           of each trial, in seconds)
            % stop   - Last time point to extract
            % step   - Step size in seconds.
            % interpolation -  enum('nearest','linear','spline','pchip','makima')
            %               Interpolation method; see timetable/synchronize. ['linear']
            % OUTPUT
            %  [t,v]  = t: time in seconds since first frame event,
            %           v: Matrix with time along the rows, and trials
            %           along columns.
            % Alternatively, when only a single output is requested:
            % T     = timetable with each column a trial. Time is in seconds
            %           relative to the first frame of the trial.
            arguments
                roi (1,1) sbx.Roi
                expt (1,1) ns.Experiment
                pv.fetchOptions = ''
                pv.modality = 'spikes'
                pv.trial = []
                pv.start =-0.5
                pv.stop  =2.5
                pv.interpolation {mustBeText} = 'linear'%, mustBeMember(pv.interpolation,{'linear','nearest','spline','pchip','makima'})}= 'linear';
                pv.step = 0.1;
            end
            if count(expt)~=1
                error("sbx.Roi.get requires a single experiment as its input.\n")
            end
            %% Get the mapping from Frames to trials.
            % Specific or this ROI (i.e. this Preprocessed set) in this
            % Expt
            trialMap = sbx.PreprocessedTrialmap & expt;
            trialMap = fetch(trialMap ,'*');
            if isempty(pv.trial)
                trials = [trialMap.trial]; % All trials
            else
                trials = pv.trial;
            end
            % Retrieve the activity in the entire session
            if ~isempty(pv.fetchOptions)
                sessionActivity= fetch(roi,pv.modality,pv.fetchOptions); % Values (e.g., spikes) across session
            else
                sessionActivity= fetch(roi,pv.modality); % Values (e.g., spikes) across session
            end
            V = [sessionActivity.(pv.modality)]; %[nrFramesPerSession nrROIs]

            newTimes = seconds(pv.start:pv.step:pv.stop);
            nrTimes  = numel(newTimes);


            nrTrials = numel(trials);
            varNames = "Trial" + string([trialMap(trials).trial]);
            T =timetable('Size',[nrTimes nrTrials],'RowTimes',newTimes','VariableTypes',repmat("doublenan",[1 nrTrials]),'VariableNames',varNames);
            trCntr=0;
            for tr = trials
                trCntr= trCntr+1;
                thisT = timetable(seconds(trialMap(tr).trialtime),V(trialMap(tr).frame,:));
                if pv.start <0 && trialMap(tr).trial>1
                    previousTrial = trialMap(tr).trial-1;
                    preTime = seconds(trialMap(previousTrial).trialtime);
                    preTime = preTime-preTime(end)-mean(diff(preTime));
                    preT = timetable(preTime,V(trialMap(previousTrial).frame,:));
                    thisT = [preT;thisT]; %#ok<AGROW>
                end
                thisT = retime(thisT,newTimes,pv.interpolation,'EndValues',NaN);
                T.(varNames(trCntr)) = table2array(thisT);
            end
            % Return as doubles or as timetable.
            if nargout ==2
                varargout{1} = seconds(T.Time);
                varargout{2} = double(T{:,1:end});
            else
                varargout{1} =T;
            end
        end
    end

    methods (Access=protected)
        function makeTuples(tbl,key)
            % Read the npy results from the suit2p folder and store them in
            % the table.
            CHUNK =350;  % This many ROIs are sent to the server at the same time.

            fldr= getFolder(sbx.Preprocessed & key);
            planes = dir(fullfile(fldr,'plane*'));
            files= {'iscell','F','Fneu','spks'};

            %% Get calibration info
            micPerPix = sqrt(sum(sbx.micronPerPixel(key).^2));

            for pl = 1:numel(planes)
                %% Read npy
                vals=  cell(1,numel(files));
                fprintf('Reading numpy files...\n')
                for f=1:numel(files)
                    thisFile = fullfile(fldr,planes(pl).name,[files{f} '.npy']);
                    if ~exist(thisFile,"file")
                        error('File %s does not exist',thisFile);
                    end
                    vals{f}= single(py.numpy.load(thisFile,allow_pickle=true));
                end
                thisFile = fullfile(fldr,planes(pl).name,'ops.npy');
                ops = py.numpy.load(thisFile,allow_pickle=true);
                % We saved the stat.npy as stat.mat in sbx.Preprocessed
                thisFile = fullfile(fldr,planes(pl).name,'stat.mat');
                load(thisFile,'stat');
                stat= [stat{:}];
                med= cat(1,stat.med); %[x y] pixels per ROI.
                compact = cat(1,stat.compact);
                aspect = cat(1,stat.aspect_ratio);

                radius = cat(1,stat.radius); % Pixels
                radius = radius.*micPerPix;


                fprintf('Done.\n')
                %% Make tuples and insert
                [iscell,f,fneu,spks] = deal(vals{:});

                [nrROIs,nrFrames] = size(f); %#ok<ASGLU>
                frameDuration = 1./ops.item{'fs'};
                meanrate = mean(spks,2,"omitnan")/frameDuration;
                stdrate = std(spks,0,2,"omitnan")/frameDuration;
                tpl = struct('subject',key.subject,...
                    'session_date',key.session_date,...
                    'prep',key.prep,...
                    'roi',num2cell(1:nrROIs)', ...
                    'plane',pl-1, ...
                    'pcell',num2cell(iscell(:,2)), ...
                    'fluorescence',num2cell(f',1)',...
                    'neuropil',num2cell(fneu',1)',...
                    'spikes',num2cell(spks',1)', ...
                    'meanrate',num2cell(meanrate',1)',...
                    'stdrate', num2cell(stdrate',1)', ...
                    'x',num2cell(med(:,1)), ...
                    'y',num2cell(med(:,2)), ...
                    'radius',num2cell(radius',1)', ...
                    'compact',num2cell(compact',1)', ...
                    'aspect',num2cell(aspect',1)' );
                fprintf('Adding %d ROIs to the database in chunks of %d\n',nrROIs,CHUNK);
                tic;
                for i=1:CHUNK:numel(tpl)
                    insert(tbl,tpl(i:min(i+CHUNK-1,numel(tpl))));
                end
                fprintf('Done in %s \n',seconds(toc))
            end
        end
    end
end