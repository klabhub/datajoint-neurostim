%{
# Stores eye position information 
-> ns.Experiment # Corresponding experiment
-> sbx.EyeParms  # The parameters that define the pose extraction process.
---
x :longblob  # The x position of the pupil center; [nrTimePoints 1]
y :longblob  # The y position of the pupil center; [nrTimePoints 1]
a: longblob  # The pupil area
quality : longblob # The quality of the estimation at each time point [nrTimePoints 1]
manualqc = NULL : smallint # quantify the overall quality based on manual inspection. 
nrtimepoints :  int unsigned # Number of time points in the pose estimation
width  : float # Width of the camera image
height : float # Height of the camera image
%}
%
% BK - Sept 2023.

classdef Eye < dj.Computed
    properties (Dependent)
        keySource
    end


    methods
        function v= get.keySource(~)
            v = (proj(ns.Experiment) & (ns.File & 'filename LIKE ''%_eye.mp4'''))*sbx.EyeParms;
        end
    end

    methods (Access=public)
        function plot(tbl,pv)
            arguments
                tbl (1,1) sbx.Eye
                pv.mode (1,1) string {mustBeMember(pv.mode,["MOVIE","TRAJECTORY","TIMECOURSE"])} = "TRAJECTORY"
            end

            for tpl = tbl.fetch('*')'
                figName= sprintf('#%s on %s@%s',tpl.subject, tpl.session_date,tpl.starttime);
                figure('Name',figName);
            switch upper(pv.mode)
                case "MOVIE"
                case "TRAJECTORY"
                    % Show trajectory x,y, area
                    scatter(tpl.x, tpl.y,tpl.a,'ko');
                    set(gca,'XLim',[1 tpl.width],'Ylim',[1 tpl.height]);
                    xlabel 'X (pixels)';
                    ylabel 'Y (pixels)';
                    case "TIMECOURSE"
                    t =1:tpl.nrtimepoints;
                    yyaxis left
                    plot(t,tpl.x./tpl.width,'LineStyle','-','Color','r');
                    hold on
                    plot(t,tpl.y/tpl.height,'LineStyle','-','Color','g');
                    set(gca,'YLim',[0 1])
                    yyaxis right;
                    plot(t,tpl.a/max(tpl.a),'LineStyle','-','Color','b');
                    hold on
                    plot(t,tpl.quality,'LineStyle','-','Color','c');
                    set(gca,'YLim',[0 1])
                    legend('x','y','area','quality')
               
            end
                       title(sprintf('#%s on %s@%s : mean quality= %.2f  NaN-Frac=%.2f',tpl.subject, tpl.session_date,tpl.starttime,mean(tpl.quality),mean(isnan(tpl.a))))
     
            end

        end
    end
    methods (Access = protected)
        function makeTuples(tbl,key)
            parms= fetch1(sbx.EyeParms &key,'parms');
            filename = fetch1(ns.File & key & 'filename LIKE ''%_eye.mp4''','filename');
            videoFile = fullfile(getenv('NS_ROOT'),filename);
            if ~exist(videoFile,'file')
                error('Eye tracker movie file %s not found',videoFile)
            end
            switch upper(key.tag)
                case 'IMFINDCIRCLES'
                    %% Pupil tracking, using imfindcircles
                    [x,y,a,quality,nrT,w,h] = sbx.Eye.imfindcircles(videoFile, parms);
                case 'DLC'

                otherwise
                    error('Unknown %d tag',key.tag);
            end
          
            tpl = mergestruct(key,struct('x',x,'y',y,'a',a,'quality',quality,'nrtimepoints',nrT,'width',w,'height',h));
            insert(tbl,tpl);
        end
    end

    methods (Static)
        function [x,y,a,quality,nrT, w,h] =imfindcircles(movieFile,pv)
            %% Use parpool if available
            if isempty(gcp('nocreate'))
                nrWorkers = 0;
            else
                nrWorkers  = gcp('nocreate').NumWorkers;
            end
            maxFrames = Inf;%30*60;  % Used for debugging only.

            % Open the movie (takes a while as it counts the frames by reading the entire movie.)
            movie = VideoReader(movieFile);     
            w = movie.Width;
            h = movie.Height;
            nrT = min(maxFrames,movie.NumFrames);

            %% Initialize output vars
            x       = nan(nrT,1);
            y       = nan(nrT,1);
            a       = nan(nrT,1);
            quality  = nan(nrT,1);

            warnNoTrace('Eye tracking analysis on %d workers\n',nrWorkers)
            % Read the movie into memory (could check that we have
            % enough...)      
            
            frames= single(movie.read([1 nrT]));
            parfor (f=1:nrT,nrWorkers)
                %for f=1:nrT  % debug
                [center,radius,thisQuality] = imfindcircles(frames(:,:,f),pv.radiusRange,'ObjectPolarity',pv.objectPolarity,'Method',pv.method,'Sensitivity',pv.sensitivity,'EdgeThreshold',pv.edgeThreshold); %#ok<PFBNS>
                if ~isempty(center)
                    [quality(f),idx] = max(thisQuality); % pick the circle with best score
                    x(f) = center(idx,1);
                    y(f) = center(idx,2);
                    a(f) = pi*radius(idx)^2;
                end
            end



        end
    end

end
