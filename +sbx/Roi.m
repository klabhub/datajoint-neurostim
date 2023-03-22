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
            %       this can be visualized.
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

            [x,y,radius,m,sd] = fetchn(roi,'x','y','radius','meanrate','stdrate');
            micPerPix = sbx.micronPerPixel(roi);  % Scaling
            mixPerPixR = sqrt(sum(micPerPix.^2));
            if isempty(pv.sz)
                % Use the physical size as the size of the cells
                pv.sz = pi*(radius*mixPerPixR).^2;
                pv.szLabel = 'size (\mum^2)';
            end

            zeroSize = pv.sz==0;
            if any(zeroSize)
                pv.sz(zeroSize) = eps;
            end


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
                [nrY,nrX] =size(meanImg);
                xImg = (0:nrX-1)*micPerPix(1);
                yImg = (0:nrY-1)*micPerPix(2);
                imagesc(xImg,yImg,meanImg)
                axImg =gca;
                axis(axImg,"equal")
                colormap(axImg,gray)
                xl = [min(xImg) max(xImg)];
                yl = [min(yImg) max(yImg)];
                axImg.XLim =xl;
                axImg.YLim = yl;
                ax = axes('position',get(axImg,'position'),'Color','none');
            else
                ax =gca;
            end

            % The data 
            hScatter = scatter(ax,y*micPerPix(1),x*micPerPix(1),pv.sz, pv.color,'filled');
            xlabel 'Y (\mum)'
            ylabel 'X (\mum)'
            % Match Suite2p layout
            set(ax,'XDir','normal','YDir','reverse');
            if ~isempty(pv.clim)
                clim(ax,pv.clim)
            end
            axis(ax,"equal")
            ax.Color = "none";
            ax.XLim =xl;
            ax.YLim = yl;
            h = colorbar;
            ylabel(h,pv.colorLabel);
            set(axImg,'position',get(ax,'position'))

            % Data tips
            if ~isempty(pv.szLabel)
                hScatter.DataTipTemplate.DataTipRows(3).Label = pv.szLabel;
            end
            if ~isempty(pv.colorLabel)
                hScatter.DataTipTemplate.DataTipRows(4).Label = pv.colorLabel;
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
                pv.modality = 'spikes'
                pv.trial = []
                pv.start =-0.5
                pv.stop  =2.5
                pv.interpolation {mustBeText} = 'linear'%, mustBeMember(pv.interpolation,{'linear','nearest','spline','pchip','makima'})}= 'linear';
                pv.step = 0.1;
            end

            %% Get the mapping from Frames to trials.
            % Specific or this ROI (i.e. this Preprocessed set) in this
            % Expt
            frame = sbx.Frame & expt;
            if ~isempty(pv.trial)
                frame = frame & struct('trial',num2cell(pv.trial)');
            end
            frames= fetch(frame,'*');

            % Retrieve the activity in the entire session
            sessionActivity= fetch(roi,pv.modality); % Values (e.g., spikes) across session
            V = [sessionActivity.(pv.modality)]; %[nrFramesPerSession nrROIs]

            % Define the new time axis (time relative to firstFrame event).
            newTimes = seconds(pv.start:pv.step:pv.stop);
            for f = frames'
                thisT = timetable(seconds(f.trialtime),V(f.frame,:),'VariableNames',"Trial" + string(f.trial));
                if f.trial ==1
                    T= retime(thisT,newTimes,pv.interpolation,'EndValues',NaN);
                else
                    T = synchronize(T,thisT, newTimes, pv.interpolation, 'EndValues',NaN);
                end
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