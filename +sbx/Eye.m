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
framerate : float # Framerate of the movie
%}
%
% BK - Sept 2023.

classdef Eye < dj.Computed
    properties (Dependent)
        keySource
    end

    methods
        function v= get.keySource(~)
            v = (proj(ns.Experiment) & (ns.File & 'filename LIKE ''%_eye.%'''))*sbx.EyeParms;
        end
    end

    methods (Access=public)
        function movie=openMovie(~,key)
            % It is possible that there is more than one _eye file; pick the smallest one
            % (presumably this is a preprocess/compressed version)
            fldr = folder(ns.Experiment& key);
            minSize = Inf;
            for f=fetch(ns.File & key & 'filename LIKE ''%_eye%''','filename')'
                    ff =fullfile(fldr,f.filename);
                    d = dir(ff);
                    if d.bytes<minSize
                        minSize= d.bytes;
                        movieFile= ff;
                    end
            end
            if ~exist(movieFile,"file")
                error('%s file not found. (Is NS_ROOT set correctly (%s)?)',movieFile,getenv('NS_ROOT'));
            end
            movie = VideoReader(movieFile);
        end

        function plot(tbl,pv)
            arguments
                tbl (1,1) sbx.Eye
                pv.mode (1,1) string {mustBeMember(pv.mode,["MOVIE","TRAJECTORY","TIMECOURSE"])} = "TRAJECTORY"
            end

            for tpl = tbl.fetch('*')'
                figName= sprintf('#%s on %s@%s',tpl.subject, tpl.session_date,tpl.starttime);
                figByName(figName);
                clf;
                switch upper(pv.mode)
                    case "MOVIE"
                        % Show the movie with the decoded pupil on top.
                        movie = openMovie(tbl,tpl);
                        frameCntr = 0;
                        phi = linspace(0,2*pi,100);
                        while (movie.hasFrame)
                            frameCntr = frameCntr+1;
                            frame = movie.readFrame;
                            hold off
                            imagesc(frame);
                            hold on
                            plot(tpl.x(frameCntr),tpl.y(frameCntr),'r*');
                            radius = sqrt(tpl.a(frameCntr)/pi);

                            line(tpl.x(frameCntr)+radius.*cos(phi),tpl.y(frameCntr)+radius.*sin(phi),'Color','g')
                            xlabel 'X (pixels)';
                            ylabel 'Y (pixels)';
                            drawnow;
                        end
                    case "TRAJECTORY"
                        % Show trajectory x,y, area on the screen.
                        scatter(tpl.x, tpl.y,tpl.a,'ko');
                        set(gca,'XLim',[1 tpl.width],'Ylim',[1 tpl.height]);
                        xlabel 'X (pixels)';
                        ylabel 'Y (pixels)';

                        title(sprintf('#%s on %s@%s : mean quality= %.2f  NaN-Frac=%.2f',tpl.subject, tpl.session_date,tpl.starttime,mean(tpl.quality,'omitnan'),mean(isnan(tpl.a))));

                    case "TIMECOURSE"
                        % Show x,y, area as a function of time.
                        T=tiledlayout(3,1,"TileSpacing","tight");
                        t =(0:tpl.nrtimepoints-1)/tpl.framerate;
                        nexttile(T)
                        plot(t,tpl.x./tpl.width);
                        hold on
                        plot(t,tpl.y/tpl.height);
                        ylim([0 1]);
                        ylabel 'Position (frac)'
                        legend('x','y')
                        nexttile(T)
                        scaledA = (tpl.a-mean(tpl.a,"omitnan"))/std(tpl.a,0,"omitnan");
                        plot(t,scaledA);
                        ylabel 'Area (z-score)'
                        legend('area')
                        nexttile(T)
                        plot(t,tpl.quality);
                        ylabel 'Position Quality ([0 1])'
                        xlabel 'Time (s)'
                        ylim([0 1]);
                        legend('quality')
                        title(T,sprintf('#%s on %s@%s : mean quality= %.2f  NaN-Frac=%.2f',tpl.subject, tpl.session_date,tpl.starttime,mean(tpl.quality,'omitnan'),mean(isnan(tpl.a))));
                end
            end
        end
    end


    methods (Access = protected)
        function makeTuples(tbl,key)
            parms= fetch1(sbx.EyeParms &key,'parms');
            movie = openMovie(tbl,key);
            switch upper(key.tag)
                case 'IMFINDCIRCLES'
                    %% Pupil tracking, using imfindcircles
                    [x,y,a,quality,nrT,w,h,fr] = sbx.Eye.imfindcircles(movie, parms);
                case 'DLC'
                    [x,y,a,quality,nrT,w,h,fr] = sbx.Eye.dlc(movie, parms);
                otherwise
                    error('Unknown %d tag',key.tag);
            end
            tpl = mergestruct(key,struct('x',x,'y',y,'a',a,'quality',quality,'nrtimepoints',nrT,'width',w,'height',h,'framerate',fr));
            insert(tbl,tpl);
        end
    end

    methods (Static)
        function [x,y,a,quality,nrT, w,h,fr] =imfindcircles(movie,pv)
            % The imfindcircles tool defines the following parameters,
            % which should be specified in the pv struct.  These are the
            % defaults (see help  imfindcircles).
            % 'objectPolarity','bright'
            % 'method','PhaseCode'
            % 'sensitivity',0.85,...
            % 'edgeThreshold',[],...
            %'radiusRange',[6 20]);
            %
            % NOTE:
            % The imfindcircles algorithm fails to find some pretty obvious
            % circles in eye tracking movies. I do not know why. Also, this
            % algorithm treats each frame on its own even though x,y, and
            % area cannot change that rapidly. It may be possible to build
            % in some filtering (e.g., search only an area around the
            % previously identified location).
            %
            arguments
                movie (1,1) VideoReader
                pv (1,1) struct
            end
            %% Use parpool if available
            if isempty(gcp('nocreate'))
                nrWorkers = 0;
            else
                nrWorkers  = gcp('nocreate').NumWorkers;
            end         

            w = movie.Width;
            h = movie.Height;
            nrT = movie.NumFrames;
            fr = movie.FrameRate;
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


        function [x,y,a,quality,nrT, w,h,fr] =dlc(movie,pv)
            % Use DeepLabCut to determine the pupil position. The actual
            % work is done by calling a DLC singularity container, and
            % parameters of dlc (e.g., which trained network) are specified
            % as pv (which is pulled from the EyeParms table in
            % makeTuples).
                

        end
    end

end
