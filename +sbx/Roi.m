%{
# Properties of a ROI in a session, based on a Preprocessed set.
-> ns.CChannel
roi : smallint    #  id within this segmentation and plane
plane : smallint #  plane
---
pcell = 0 : float # Probability that the ROI is a neuron
x        : Decimal(4,0) # x-Pixel location 
y        : Decimal(4,0) # y-Pixel location 
radius   : Decimal(4,1) # Radius of cell in micron
aspect   : Decimal(4,1)  # Aspect raio of major/minor axis of a 2D Gaussian fit.
compact  : Decimal(4,2)  # How compact the ROI is ( 1 is a disk, >1 means less compact)
%}
classdef Roi < dj.Computed
    properties (Dependent)
        keySource
    end
    methods
        function v = get.keySource(~)
            % The key source is the ns.C table becuase we read all channels
            % for a ns.C entry from a single npy file
            v = ns.C & (ns.CParm & 'extension=''.sbx''');
        end
    end
    methods  (Access=public)
        function [hData] = plotSpatial(roi,pv)
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
            % EXAMPLE
            %
            %
            % plot(sbx.Roi); % Encode all cell radius (size) and z-scored session activity (color).
            % z = computed properties of all rois
            % plot(sbx.Roi,color =z,clim = [0 5]);  % Encode z by the color, crop colors at 5
            arguments
                roi (1,1) sbx.Roi
                pv.sz   (1,:)  double  = []
                pv.szLabel (1,1) string = ''
                pv.color   (1,:) double = []
                pv.colorLabel (1,1) string = ''
                pv.clim    (1,:) double  =[];
                pv.showImg (1,1) logical = true;
                pv.colormap (:,3) double = hot;
                pv.pix (1,1) logical  = false
                pv.alpha (1,1) double = 0.5;
                pv.alphaThreshold (1,1) double =0;
            end

            if count(roi) ==0
                fprintf('No rois, nothing to plot \n');
                return;
            end

            % To visualize an arbitrary and large subset of rois, this
            % function can be called with a large struct array (one per roi),
            % which is potentially inefficient. Silencing DJ warning about this.
            warnState = warning('query');
            warning('off', 'DataJoint:longCondition');

            MAXPOINTS = 300; %
            [x,y,radius,roiNr] = fetchn(roi,'x','y','radius','roi');
            [xScale,yScale] =fetchn(sbx.Preprocessed & roi,'xscale','yscale');
            micPerPix = [xScale yScale];
            micPerPixR = sqrt(sum(micPerPix.^2));


            %% Setup color
            if isempty(pv.color)
                % Use the z-scored rate as the color of the cells
                z = m./sd;
                pv.color = z;
                pv.colorLabel = 'rate (Z)';
            end
            % Clamp to the limits
            if isempty(pv.clim)
                pv.clim = [min(pv.color) max(pv.color)];
            else
                pv.color = max(min(pv.color,pv.clim(2)),pv.clim(1)); % Clamp between clims
            end
            colormap(pv.colormap)


            %% Show mean image as background
            if pv.showImg
                meanImg = fetch1(sbx.Preprocessed & roi,'img','LIMIT 1');
                RI = imref2d(size(meanImg),micPerPix(1),micPerPix(2));
                cl = [min(meanImg,[],"all") max(meanImg,[],"all")];
                meanImg  = 1+round((size(gray,1)-1)*(meanImg-cl(1))./(cl(2)-cl(1)));
                % Convert to true color to reuse axes for color map
                meanImg = ind2rgb(meanImg,gray);
                imshow(meanImg,RI);
                hold on
            end

            if pv.pix
                %% Show the data in the ROI pixels
                overlay = zeros(RI.ImageSize);
                alpha   = zeros(RI.ImageSize);
                prep = sbx.Preprocessed & roi;
                stat =  prep.stat; % Reads the stat.npy file
                assert(~isempty(stat),'Stat file not loaded. Cannot use pix=true');
                
                rCntr  =0;
                % Scale to the full colormap
                scaled  = (pv.color-min(pv.color))./(max(pv.color)-min(pv.color));
                clrIndex  = 1+round((size(pv.colormap,1)-1)*scaled);
                for r=roiNr'
                    rCntr= rCntr+1;
                    x =stat(r).xpix;
                    y = stat(r).ypix;
                    ix =sub2ind(size(overlay),y,x);
                    overlay(ix)= clrIndex(rCntr);
                    alpha(ix) = pv.alpha*(abs(pv.color(rCntr))>pv.alphaThreshold);
                end
                overlay = ind2rgb(overlay,pv.colormap); % Convert to true color using the colormap
                hData = imshow(overlay,RI);
                set(hData,'AlphaData',alpha)
            else
                %% Show the data as circles
                % Setup size
                if isempty(pv.sz)
                    % Use the physical size as the size of the cells
                    pv.sz =min(pi*(radius*micPerPixR).^2,MAXPOINTS);
                    pv.szLabel = 'size (\mum^2)';
                end
                zeroSize = pv.sz==0;
                if any(zeroSize)
                    pv.sz(zeroSize) = eps;
                end

                % Show the data symbolically on top of the image. Match
                % orientation to suite2p
                hData = scatter((x-1)*micPerPix(1),(y-1)*micPerPix(2),pv.sz(:), pv.color(:),'filled','MarkerFaceAlpha',pv.alpha);%.*(abs(pv.color(:))>pv.alphaThreshold));
                set(gca,'Color','none','XDir','normal','YDir','reverse','colormap',pv.colormap);

                % Data tips
                if ~isempty(pv.szLabel)
                    hData.DataTipTemplate.DataTipRows(3).Label = pv.szLabel;
                end
                if ~isempty(pv.colorLabel)
                    hData.DataTipTemplate.DataTipRows(4).Label = pv.colorLabel;
                end
            end
            if ~isempty(pv.clim)
                clim(pv.clim)
            end

            %% Add Markup
            xlabel 'X (\mum)'
            ylabel 'Y (\mum)'
            h = colorbar;
            ylabel(h,pv.colorLabel);


            % Restore warning state
            warning(warnState);

        end
    end
    methods (Access=protected)
        function makeTuples(tbl,cKey)
            % Extract info from Channel key
            c = fetch((ns.C & cKey)*(ns.CParm & 'extension=''.sbx'''),'*');
            if isempty(c)
                error('This ns.C (%s) does not contain SBX data. Cannot create sbx.ROI tuples.',cKey.ctag)
            end
            prep = sbx.Preprocessed & struct('subject',c.subject,'session_date',c.session_date,'prep',c.parms.prep);
            info = sbx.readInfoFile(fetch(ns.Experiment&cKey,'LIMIT 1'));
            micPerPix = sqrt(sum([info.xscale info.yscale].^2));
            if count(prep)~=1
                error('Need exactly 1 sbx.Preprocessed set to extract ROIs, not %d',count(prep))
            end
            fldr = fullfile(folder(ns.Experiment & cKey),fetch1(prep,'folder'));
            if ~exist(fldr,"dir")
                error('Preprocessed data folder %s not found',fldr)
            end
            planes = dir(fullfile(fldr,'plane*'));
            channelsSoFar = 0;

            for pl = 1:numel(planes)
                %% Read npy
                thisFile = fullfile(fldr,planes(pl).name,'iscell.npy');
                if ~exist(thisFile,"file")
                    error('File %s does not exist',thisFile);
                end
                iscell = single(py.numpy.load(thisFile,allow_pickle=true));

                % We saved the stat.npy as stat.mat in sbx.Preprocessed
                thisFile = fullfile(fldr,planes(pl).name,'stat.mat');
                if ~exist(thisFile,"file")
                    error('File %s does not exist',thisFile);
                end
                load(thisFile,'stat');
                stat= [stat{:}];
                med= cat(1,stat.med); %[y x] center pixels per ROI.
                compact = cat(1,stat.compact);
                aspect = cat(1,stat.aspect_ratio);

                radius = cat(1,stat.radius); % Pixels
                radius = radius.*micPerPix;

                %% Make tuples and insert=
                nrROIs= numel(aspect);
                key = repmat(ns.stripToPrimary(ns.C,cKey),[nrROIs 1]);
                channelsPerRoi = num2cell(channelsSoFar+(1:nrROIs));
                [key.channel] = deal(channelsPerRoi{:});
                tpl = mergestruct(key, ...
                    struct('roi',num2cell(1:nrROIs)', ...
                    'plane',pl-1, ...
                    'pcell',num2cell(iscell(:,2)), ...
                    'x',num2cell(med(:,2)), ...    % Dim 2 is the horizontal axis of the image
                    'y',num2cell(med(:,1)), ...    % Dim 1 is the vertical axis of the image
                    'radius',num2cell(radius',1)', ...
                    'compact',num2cell(compact',1)', ...
                    'aspect',num2cell(aspect',1)'));
                insert(tbl,tpl);
                channelsSoFar = channelsSoFar + nrROIs; % Next plane starts at nrROIs+1
            end
        end
    end
end