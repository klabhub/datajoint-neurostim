%{
# Stores ball velocityinformation 
-> ns.Movie     # Corresponding movie
-> sbx.BallParms  # The parameters that define the extraction process.
---
velocity :longblob  # The velocity; complex number representing the instantanous velocity of the mouse.[nrTimePoints 2]
quality : longblob # The quality of the estimation at each time point [nrTimePoints 1]
manualqc = NULL : smallint # quantify the overall quality based on manual inspection. 
%}
%
% For an experiment with two _ball files, the makeTuples function will use the
% smalllest of the two. This handles the case where large files from SBX
% were compressed after the experiment into a more tractable format/size.
%
% If a GPU is available this code will use it to speed up xcorr2.
%
% The parms struct (set in sbx.BallParms) needs two parameters:
% .scaleFactor  By how much to scale the raw images (0.5 means reduce the size)
% .minPixels    The minimum pixels along either width or height after
%               scalng by scaleFactor (in other words if scaleFactor is too small, it
%               will be increased such that minPixels remain).
%
%
% BK - Sept 2023.

classdef Ball < dj.Computed
    properties (Dependent)
        keySource
    end


    methods
        function v= get.keySource(~)
            v = (ns.Movie & 'filename LIKE ''%_ball.mj2''')*sbx.BallParms;
        end

    end

    methods (Access=public)
        function plot(tbl,pv)
            % Function to plot sbx.Ball movies for each row in the tbl.
            arguments
                tbl (1,1) sbx.Ball
                pv.mode (1,1) string {mustBeMember(pv.mode,["MOVIE","TRAJECTORY","TIMECOURSE"])} = "TRAJECTORY"
                pv.history (1,1) double = 5         % How many vectors to show recent directions
                pv.frameStep (1,1) double {mustBeInteger,mustBeNonnegative} =1
                pv.frameStart (1,1) double {mustBeInteger,mustBeNonnegative} =1
                pv.frameStop (1,1) double {mustBeInteger,mustBeNonnegative} = 100000000
            end
            tbl = tbl*ns.Movie;
            for tpl = tbl.fetch('*')'
                figName= sprintf('#%s on %s@%s',tpl.subject, tpl.session_date,tpl.starttime);
                hFig = figByName(figName);
                hFig.Units= 'Normalized';
                clf;
                switch upper(pv.mode)
                    case "MOVIE"
                        fprintf('Opening movie file...\n')
                        movie = open(ns.Movie & tpl,"smallest");
                        fprintf('done.\n Press Ctrl-C to stop, or close the window to move to the next movie in the table\n');
                        ax = axes('Position',[0 0 1 1]);
                        axis(ax,'off')
                        pos = get(ax,'Position');
                        axPolar = polaraxes(hFig,'Position',[0 0 pos(3:4)/5]);
                        axSpeed = axes(hFig,'Position',[0.75 0.05 pos(3:4)/5]);
                        xlabel(axSpeed,'Frame #','Color','y','FontWeight','Bold','FontSize',12);
                        ylabel(axSpeed,'Speed','Color','y','FontWeight','Bold','FontSize',12)
                        speed =abs(tpl.velocity);
                        maxSpeed = max(eps,max(speed));
                        historyColormap = gray; % Show multiple trailing vectors as shades of grays
                        colormap gray
                        for frameCntr = pv.frameStart:pv.frameStep:min(pv.frameStop,tpl.nrframes)
                            if  ~ishandle(hFig);break;end;
                            frame = movie.read(frameCntr);
                            hold off
                            imagesc(ax,frame);
                            hold on
                            text(ax,max(ax.XLim), min(ax.YLim),sprintf('%s       %dx - Frame #%d/%d',movie.Name,pv.frameStep,frameCntr,tpl.nrframes),'HorizontalAlignment','Right','VerticalAlignment','top','Color','y','FontWeight','Bold','FontSize',12,'Interpreter','none')
                            % Instantaneous velocity with trailing vectors
                            % for history
                            if ~isnan(tpl.velocity(frameCntr))
                                fToKeep = frameCntr-pv.history:frameCntr;
                                fToKeep(fToKeep<1) =[];
                                nrF =numel(fToKeep);
                                theta  = angle(tpl.velocity(fToKeep));
                                rho   = speed(fToKeep);
                                h = polarplot(axPolar,[zeros(1,nrF);theta'],[zeros(1,nrF);rho']);
                                % Use shading such that the most recent frame
                                % is black and earlier ones fade to white.
                                ix = round(linspace(255,1,nrF));

                                colors = num2cell(historyColormap(ix,:),2)';
                                [h.Color] =deal(colors{:});
                                axPolar.RLim = [0 1.1];
                            end

                            %% speed with 100 frames history
                            nrFramesToKeep= 100;
                            fToKeep = frameCntr+(-nrFramesToKeep:0);
                            fToKeep(fToKeep<1)=[];
                            plot(axSpeed,fToKeep,speed(fToKeep),'r');
                            xlim(axSpeed,[frameCntr-nrFramesToKeep frameCntr])
                            set(axSpeed,'YTick',[],'YScale','Linear')
                            ylim(axSpeed,[eps 1.1*maxSpeed]);
                            axSpeed.XColor ='y';


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
                        % Show trajectory by "integrating" over time (cumsum).
                        % Starting point is always 0,0
                        % Time is shown as hsv color
                        % Speed is shown as marker size

                        stay =~any(isnan(tpl.velocity),2);
                        position  = cumsum(tpl.velocity(stay,:));
                        speed  = 1+abs(tpl.velocity(stay,:));
                        time= hsv(sum(stay)); % Show multiple trailing vectors as shades of grays
                        scatter(0,0,10,'*');
                        hold on
                        scatter(real(position),imag(position),speed,time,'o','filled');
                        xlabel 'X (pixels)';
                        ylabel 'Y (pixels)';
                        set(gca,'XLim',max(abs(xlim))*[-1 1],'ylim',max(abs(ylim))*[-1 1]);
                        title(sprintf('#%s on %s@%s : mean quality= %.2f  NaN-Frac=%.2f',tpl.subject, tpl.session_date,tpl.starttime,mean(tpl.quality,'omitnan'),mean(isnan(tpl.velocity))));
                        colormap hsv
                        h = colorbar;
                        ylabel(h,'Time (norm)')
                        set(h,'YTick',0:0.25:1)
                    case "TIMECOURSE"
                        % Show dx,dy and quality over time.
                        T=tiledlayout(2,1,"TileSpacing","tight");
                        t =(0:tpl.nrframes-1)/tpl.framerate;
                        nexttile(T)
                        plot(t,real(tpl.velocity));
                        hold on
                        plot(t,imag(tpl.velocity));
                        ylabel 'Speed (pixels/frame)'
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
            mvFile = file(ns.Movie &key,"smallest");
            [pth,filename] = fileparts(mvFile);
            csvFile =fullfile(pth,filename + '_' + key.tag + '.csv');
            if exist(csvFile,"file")
                fprintf('%s already exists. Adding to the table.\n',csvFile);
                T=readtable(csvFile,'ReadVariableNames',true);
                velocity = T.x +1i.*T.y;
                tpl = mergestruct(key,struct('velocity',velocity,'quality',T.quality));
            else
                movie = open(ns.Movie &key,"smallest");
                switch upper(key.tag)
                    case 'XCORR'
                        [velocity,quality] =sbx.Ball.xcorr(movie,parms);
                    case 'PHASECORR'
                        [velocity,quality] =sbx.Ball.phasecorr(movie,parms);
                    otherwise
                        error('Unknown %d tag',key.tag);
                end
                tpl = mergestruct(key,struct('velocity',velocity,'quality',quality));

                %%  Also save the results in a csv file
                T = table((1:numel(quality))',real(velocity),imag(velocity),quality,'VariableNames',{'frame','x','y','quality'});
                writetable(T,csvFile,Delimiter=',',LineEnding='\r\n',WriteVariableNames=true,WriteMode='overwrite');
            end

            if count (tbl&ns.stripToPrimary(tbl,tpl))==0
                insert(tbl,tpl);
            else
                fprintf('%s/%s already exists in the table\n',key.filename,key.tag)
            end
        end
    end
    methods (Static)

        function [velocity,quality] =phasecorr(movie,parms)
            arguments
                movie (1,1) VideoReader
                parms (1,1) struct
            end
            useGPU = false; %canUseGPU;  % If a GPU is available we'll use it.

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

            warnNoTrace('Ball tracking phasecorr analysis (useGPU: %d)\n',useGPU);
            f=1;
            z1 = im2single(movie.readFrame);
            if ndims(z1)==3;z1=z1(:,:,1);end  % Images are gray scale but some have been saved with 3 planes of identical bits.

            % Determine the scaling.
            heightWidth  = round(parms.scaleFactor.*[h w]);
            if any(heightWidth<parms.minPixels)
                % Too small: correct
                parms.scaleFactor = max(parms.minPixels./[h w]);
                heightWidth = round(parms.scaleFactor*[w h]);
            end
            z1 =imresize(z1,heightWidth);
            while movie.hasFrame
                z2 =  im2single(movie.readFrame);
                if ndims(z2)==3;z2=z2(:,:,1);end
                z2 =imresize(z2,heightWidth);
                if useGPU
                    z1 = gpuArray(z1);
                    z2 = gpuArray(z2);
                end

                [tform,quality(f)] = imregcorr(z2,z1,"translation", "window",true);

                dy = tform.Translation(2);
                dx = tform.Translation(1);

                velocity(f) = dx  - 1i.*dy; % Reflect the motion of the mouse,not the ball
                % next frame
                f=f+1;
                z1=z2;
            end
            if useGPU
                velocity = gather(velocity);
                quality = gather(quality);
            end

        end

        function [velocity,quality] =xcorr(movie,parms)
            arguments
                movie (1,1) VideoReader
                parms (1,1) struct
            end
            useGPU = canUseGPU;  % If a GPU is available we'll use it.

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
            z1 = single(movie.readFrame);
            if ndims(z1)==3;z1=z1(:,:,1);end  % Images are gray scale but some have been saved with 3 planes of identical bits.

            % Determine the scaling.
            heightWidth  = round(parms.scaleFactor.*[h w]);
            if any(heightWidth<parms.minPixels)
                % Too small: correct
                parms.scaleFactor = max(parms.minPixels./[h w]);
                heightWidth = round(parms.scaleFactor*[w h]);
            end
            z1 =imresize(z1,heightWidth);
            while movie.hasFrame
                z2 =  single(movie.readFrame);
                if ndims(z2)==3;z2=z2(:,:,1);end
                z2 =imresize(z2,heightWidth);
                if useGPU
                    z1 = gpuArray(z1);
                    z2 = gpuArray(z2);
                end
                %% Find maximum xcorr
                %  Must substract the mean.
                z1 = z1-mean(z1,"all","omitnan");
                z2 = z2-mean(z2,"all","omitnan");
                xc =xcorr2(z1,z2);
                [maxXC,ix] = max(xc(:));
                scale=sum(((z1+z2)/2).^2,'all');
                quality(f) = maxXC./scale;
                [dy,dx]= ind2sub(size(xc),ix);
                dy = dy-heightWidth(1);
                dx = dx-heightWidth(2);
                velocity(f) = dx  - 1i.*dy; % Reflect the motion of the mouse,not the ball
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