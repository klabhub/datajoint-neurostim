%{
# Stores ball velocityinformation 
-> ns.Experiment # Corresponding experiment
-> sbx.BallParms  # The parameters that define the extraction process.
---
velocity :longblob  # The velocity; first column is the horizontal component, second column the vertical component. [nrTimePoints 2]
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
            v = (proj(ns.Experiment) & (ns.File & 'filename LIKE ''%_ball.avi'''))*sbx.BallParms;
        end

    end

    methods (Access=public)
        function movie=openMovie(~,key)
            filename = fetch1(ns.File & key & 'filename LIKE ''%_ball.avi''','filename');
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
                pv.history (1,1) double = 10                
            end

            for tpl = tbl.fetch('*')'
                figName= sprintf('#%s on %s@%s',tpl.subject, tpl.session_date,tpl.starttime);
                hFig = figByName(figName);
                hFig.Units= 'Normalized';
                clf;
                switch upper(pv.mode)
                    case "MOVIE"
                        movie = openMovie(tbl,tpl);
                        ax = axes('Position',[0 0 1 1]);
                        axis(ax,'off')
                        pos = get(ax,'Position');
                        axPolar = polaraxes(hFig,'Position',[0 0 pos(3:4)/5]);
                        axSpeed = axes(hFig,'Position',[0.75 0.05 pos(3:4)/5]);
                        maxVelocity = max(eps,prctile(abs(tpl.velocity),95));
                        historyColormap = gray; % Show multiple trailing vectors as shades of grays
                        for frameCntr = 1:tpl.nrtimepoints                            
                            frame = movie.readFrame;
                            hold off
                            imagesc(ax,frame);                             
                            hold on
                            text(ax,max(ax.XLim), min(ax.YLim),sprintf('Frame #%d/%d',frameCntr,tpl.nrtimepoints),'HorizontalAlignment','Right','VerticalAlignment','top','Color','y','FontWeight','Bold','FontSize',12)
                            % Instantaneous velocity with trailing vectors
                            % for history
                            if ~isnan(tpl.velocity(frameCntr))
                                fToKeep = frameCntr-pv.history:frameCntr;
                                fToKeep(fToKeep<1) =[];
                                nrF =numel(fToKeep);
                                h = polarplot(axPolar,[complex(zeros(1,nrF)) ;tpl.velocity(fToKeep)']);
                                % Use shading such that the most recent frame
                                % is black and earlier ones fade to white.
                                ix = round(linspace(255,1,nrF));
                                
                                colors = num2cell(historyColormap(ix,:),2)';
                                [h.Color] =deal(colors{:});                               
                                axPolar.RLim = [0 maxVelocity];
                            end

                            %% speed with 100 frames history
                            nrFramesToKeep= 100;
                            fToKeep = frameCntr+(-nrFramesToKeep:0);
                            fToKeep(fToKeep<1)=[];
                            plot(axSpeed,fToKeep,abs(tpl.velocity(fToKeep)),'r');
                            xlim(axSpeed,[frameCntr-nrFramesToKeep frameCntr])
                            ylim(axSpeed,[0 1.1*maxVelocity]);
                            set(axSpeed,'YTick',[])


                            %% Quality indicated by the circle in the polar plot. (red = bad, green is good)
                            if ~isnan(tpl.quality(frameCntr))
                                hold(axPolar,"on")
                                polarplot(axPolar,0,0,'.','MarkerSize',10,'Color',[max(0,1-tpl.quality(frameCntr)) min(1,tpl.quality(frameCntr)) 0])
                                hold(axPolar,"off")
                            else
                                polarplot(axPolar,0,0,'X','MarkerSize',10,'Color','r')
                            end
                            % polarplot puts ticks back every time.
                            axPolar.RTickLabel = [];
                            axPolar.ThetaTickLabel  = [];
                            axPolar.RTick = [];
                            axPolar.ThetaTick  = [];

                            drawnow;
                        end

                    case "TRAJECTORY"
                        % Show trajectory x,y, area
                        stay =~any(isnan(tpl.velocity),2);
                        position  = cumsum(tpl.velocity(stay,:));
                        speed  = 1+abs(tpl.velocity(stay,:));
                        time= hsv(sum(stay)); % Show multiple trailing vectors as shades of grays
                        
                        scatter(real(position),imag(position),speed,time,'o','filled');
                        xlabel 'X (a.u.)';
                        ylabel 'Y (a.u.)';
                        title(sprintf('#%s on %s@%s : mean quality= %.2f  NaN-Frac=%.2f',tpl.subject, tpl.session_date,tpl.starttime,mean(tpl.quality,'omitnan'),mean(isnan(tpl.velocity))));

                    case "TIMECOURSE"
                        T=tiledlayout(2,1,"TileSpacing","tight");
                        t =(0:tpl.nrtimepoints-1)/tpl.framerate;
                        nexttile(T)
                        plot(t,real(tpl.velocity));
                        hold on
                        plot(t,imag(tpl.velocity));                        
                        ylabel 'Speed (a.u.)'
                        legend('dx','dy')
                        nexttile(T)
                        plot(t,tpl.quality);
                        ylabel 'Velocity Quality ([0 1])'
                        xlabel 'Time (s)'
                        ylim([0 1.1]);
                        legend('quality')
                        title(T,sprintf('#%s on %s@%s : mean quality= %.2f  NaN-Frac=%.2f',tpl.subject, tpl.session_date,tpl.starttime,mean(tpl.quality,'omitnan'),mean(isnan(tpl.velocity))));
                        linkaxes(T.Children,'x')
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
        function [velocity,quality,nrFrames,fr] =xcorr(movie,pv)
            arguments
                movie (1,1) VideoReader
                pv (1,1) struct
            end
            maxFrames= Inf;%Debugging 30*60;

            useGPU = canUseGPU;

            w=movie.Width;h =movie.Height;
            nrFrames = movie.NumFrames;
            fr = movie.FrameRate;
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
            pixelScaleDown = 3;
            z1 = single(movie.readFrame);
            if ndims(z1)==3;z1=z1(:,:,1);end  % Images are gray scale but some have been saved with 3 planes of identical bits.
            z1 =imresize(z1,round([w h]./pixelScaleDown));
            [scaledH,scaledW] = size(z1); % After scaling
            while movie.hasFrame &&  f<maxFrames
                z2 =  single(movie.readFrame);
                if ndims(z2)==3;z2=z2(:,:,1);end % Force gray scale
                z2 =imresize(z2,round([w h]./pixelScaleDown));
                if useGPU
                    z1 = gpuArray(z1);
                    z2 = gpuArray(z2);
                end
                %% Find maximum xcorr
                xc =xcorr2(z1,z2);
                [maxXC,ix] = max(xc(:));
                scale=sum(((z1+z2)/2).^2,'all');
                quality(f) = maxXC./scale;
                [dy,dx]= ind2sub(size(xc),ix);
                dy = dy-scaledH;
                dx = dx-scaledW;
                velocity(f) = -dx  - 1i.*dy; % Minus sign to reflect the motion of the mouse, which is opposite to that of the ball.
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
