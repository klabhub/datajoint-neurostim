%{
# Stores ball velocityinformation 
-> ns.Experiment # Corresponding experiment
-> sbx.BallParms  # The parameters that define the extraction process.
---
v :longblob  # The velocity; first column is the horizontal component, second column the vertical component. [nrTimePoints 2]
quality : longblob # The quality of the estimation at each time point [nrTimePoints 1]
manualqc = NULL : smallint # quantify the overall quality based on manual inspection. 
nrtimepoints :  int unsigned # Number of time points in the pose estimation
framerate : float # Framerate of the movie
%}
%
% BK - Sept 2023.

classdef Ball < dj.Computed
    properties (Dependent)
        keySource
    end


    methods
        function v= get.keySource(~)
            v = (proj(ns.Experiment) & (ns.File & 'filename LIKE ''%_ball.mp4'''))*sbx.BallParms;
        end

    end

    methods (Access=public)
        function movie=openMovie(~,key)
            filename = fetch1(ns.File & key & 'filename LIKE ''%_ball.mp4''','filename');
            movieFile = fullfile(getenv('NS_ROOT'),filename);
            if ~exist(movieFile,"file")
                error('%s file not found. (Is NS_ROOT set correctly (%s)?)',movieFile,getenv('NS_ROOT'));
            end
            movie = VideoReader(movieFile);
        end

        function plot(tbl,pv)
            arguments
                tbl (1,1) sbx.Ball
                pv.mode (1,1) string {mustBeMember(pv.mode,["MOVIE","TRAJECTORY","TIMECOURSE"])} = "TRAJECTORY"
            end

            for tpl = tbl.fetch('*')'
                figName= sprintf('#%s on %s@%s',tpl.subject, tpl.session_date,tpl.starttime);
                figByName(figName);
                clf;
                switch upper(pv.mode)
                    case "MOVIE"
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
                        % Show trajectory x,y, area
                        scatter(tpl.x, tpl.y,tpl.a,'ko');
                        set(gca,'XLim',[1 tpl.width],'Ylim',[1 tpl.height]);
                        xlabel 'X (pixels)';
                        ylabel 'Y (pixels)';

                        title(sprintf('#%s on %s@%s : mean quality= %.2f  NaN-Frac=%.2f',tpl.subject, tpl.session_date,tpl.starttime,mean(tpl.quality,'omitnan'),mean(isnan(tpl.a))));

                    case "TIMECOURSE"
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
            parms= fetch1(sbx.BallParms &key,'parms');
            movie = openMovie(tbl,key);
            switch upper(key.tag)
                case 'XCORR'
                    [velocity,quality,nrT,fr] =sbx.Ball.xcorr(movie,parms);
                otherwise
                    error('Unknown %d tag',key.tag);
            end

            tpl = mergestruct(key,struct('velocity',velocity,'quality',quality,'nrtimepoints',nrT,'framerate',fr));
            insert(tbl,tpl);
        end
    end

    methods (Static)
        function [velocity,quality,nrT,fr] =xcorr(movie,pv)
            arguments
                movie (1,1) VideoReader
                pv (1,1) struct
            end
            maxFrames= 30*60;

            useGPU = canUseGPU;

            w=movie.Width;h =movie.Height;
            nrFrames = movie.NumFrames;
            %% Initialize output vars
            if useGPU
                velocity = nan(nrFrames,1,"gpuArray");
                quality  = nan(nrFrames,1,"gpuArray");
            else
                velocity = nan(nrFrames,1);
                quality  = nan(nrFrames,1);
            end

            warnNoTrace('Ball tracking analysis (useGPU: %d)\n',useGPU);
            f=1;
            z1 = single(movie.readFrame);
            if ndims(z1)==3;z1=z1(:,:,1);end
            while movie.hasFrame &&  f<maxFrames 
                z2 = single(movie.readFrame);
                if ndims(z2)==3;z2=z2(:,:,1);end
                if useGPU
                    z1 = gpuArray(z1);
                    z2 = gpuArray(z2);
                end
                
                xc =xcorr2(z1,z2);
                [maxXC,ix] = max(xc(:));
                scale=sum(((z1+z2)/2).^2,'all');
                quality(f) = maxXC./scale;
                [dy,dx]= ind2sub(size(xc),ix);
                dy = dy-h;
                dx = dx-w;
                velocity(f) = dx  + 1i.*dy;
                % next frame
                f=f+1;
                z1=z2;                
            end
            if useGPU
                velocity = gather(velocity);
                quality = gather(quality);
            end


        end
    end

end
