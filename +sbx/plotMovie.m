function plotMovie(tblOrKeys,pv)
% Show SBX _eye and _ball movies
% BK - Sept 2023.
arguments
    tblOrKeys
    pv.ball (1,1) logical = true  % Set to false to show the eye movie

    pv.mode (1,1) string {mustBeMember(pv.mode,["MOVIE","TRAJECTORY","TIMECOURSE"])} = "TRAJECTORY"
    pv.history (1,1) double = 5         % How many vectors to show recent directions
    pv.frameStep (1,1) double {mustBeInteger,mustBeNonnegative} =1
    pv.frameStart (1,1) double {mustBeInteger,mustBeNonnegative} =1
    pv.frameStop (1,1) double {mustBeInteger,mustBeNonnegative} = 100000000
end

if pv.ball
    include = '%_ball.mj2';
else
    include = '%eye.mj2';
end

% Create a table of experiment kes
tbl = ns.Experiment & proj(tblOrKeys);
for exptTpl = tbl.fetch('*')'
    fldr = folder(ns.Experiment &exptTpl);
    exptName= sprintf('#%s on %s@%s',exptTpl.subject, exptTpl.session_date,exptTpl.starttime);
    %% Check that we have preprocessed data for this Experiment
    thisC= ns.C & exptTpl & ['filename LIKE ''' include ''''];
    if exists(thisC)
        cTpl = fetch(thisC,'*');  % The entry in the C table
        cChannelTpl = fetch(ns.CChannel & ns.stripToPrimary(ns.CChannel,cTpl),'*'); % Struct array with continuous data, one struct per channel
        hFig = figByName(exptName);
        hFig.Units= 'Normalized';
        clf;
        if pv.ball
            %% BALL
            switch upper(pv.mode)
                case "MOVIE"
                    fprintf('Opening movie file...\n')
                    movie = VideoReader(fullfile(fldr,cTpl.filename)); %#ok<TNMLP>
                    fprintf('done.\n Press Ctrl-C to stop, or close the window to move to the next movie in the table\n');
                    ax = axes('Position',[0 0 1 1]);
                    axis(ax,'off')
                    pos = get(ax,'Position');
                    axPolar = polaraxes(hFig,'Position',[0 0 pos(3:4)/5]);
                    axSpeed = axes(hFig,'Position',[0.75 0.05 pos(3:4)/5]);
                    xlabel(axSpeed,'Frame #','Color','y','FontWeight','Bold','FontSize',12);
                    ylabel(axSpeed,'Speed','Color','y','FontWeight','Bold','FontSize',12)
                    speed =abs(cChannelTpl.velocity);
                    maxSpeed = max(eps,max(speed));
                    historyColormap = gray; % Show multiple trailing vectors as shades of grays
                    colormap gray
                    for frameCntr = pv.frameStart:pv.frameStep:min(pv.frameStop,cTpl.nrsamples)
                        if  ~ishandle(hFig);break;end
                        frame = movie.read(frameCntr);
                        hold off
                        imagesc(ax,frame);
                        hold on
                        text(ax,max(ax.XLim), min(ax.YLim),sprintf('%s       %dx - Frame #%d/%d',movie.Name,pv.frameStep,frameCntr,cTpl.nrsamples),'HorizontalAlignment','Right','VerticalAlignment','top','Color','y','FontWeight','Bold','FontSize',12,'Interpreter','none')
                        % Instantaneous velocity with trailing vectors
                        % for history
                        if ~isnan(velocity(frameCntr))
                            fToKeep = frameCntr-pv.history:frameCntr;
                            fToKeep(fToKeep<1) =[];
                            nrF =numel(fToKeep);
                            theta  = angle(cChannelTpl.velocity(fToKeep));
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
                        if ~isnan(cChannelTpl.quality(frameCntr))
                            hold(axPolar,"on")
                            polarplot(axPolar,0,0,'.','MarkerSize',10,'Color',[max(0,1-cChannelTpl.quality(frameCntr)) min(1,cChannelTpl.quality(frameCntr)) 0])
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

                    stay =~any(isnan(cChannelTpl.velocity),2);
                    position  = cumsum(cChannelTpl.velocity(stay,:));
                    speed  = 1+abs(cChannelTpl.velocity(stay,:));
                    time= hsv(sum(stay)); % Show multiple trailing vectors as shades of grays
                    scatter(0,0,10,'*');
                    hold on
                    scatter(real(position),imag(position),2,time,'o','filled');
                    xlabel 'X (pixels)';
                    ylabel 'Y (pixels)';
                    set(gca,'XLim',max(abs(xlim))*[-1 1],'ylim',max(abs(ylim))*[-1 1]);
                    title(sprintf('#%s on %s@%s : mean quality= %.2f  NaN-Frac=%.2f',exptTpl.subject, exptTpl.session_date,exptTpl.starttime,mean(cChannelTpl.quality,'omitnan'),mean(isnan(cChannelTpl.velocity))));
                    colormap hsv
                    h = colorbar;
                    ylabel(h,'Time (norm)')
                    set(h,'YTick',0:0.25:1)
                case "TIMECOURSE"
                    % Show dx,dy and quality over time.
                    T=tiledlayout(2,1,"TileSpacing","tight");
                    t =linspace(cTpl.time(1),cTpl.time(end),cTpl.time(3));
                    nexttile(T)
                    plot(t,real(cChannelTpl.velocity));
                    hold on
                    plot(t,imag(cChannelTpl.velocity));
                    ylabel 'Speed (pixels/frame)'
                    legend('dx','dy')
                    nexttile(T)
                    plot(t,exptTpl.quality);
                    ylabel 'Velocity Quality ([0 1])'
                    xlabel 'Time (s)'
                    ylim([0 1.1]);
                    legend('quality')
                    title(T,sprintf('#%s on %s@%s : mean quality= %.2f  NaN-Frac=%.2f',exptTpl.subject, exptTpl.session_date,exptTpl.starttime,mean(cChannelTpl.quality,'omitnan'),mean(isnan(cChannelTpl.velocity))));
                    linkaxes(T.Children,'x')

            end
        else
            %% EYE
            q = channelByName(cChannelTpl,'quality');
            x =channelByName(cChannelTpl,'x');
            y = channelByName(cChannelTpl,'y');
            a = channelByName(cChannelTpl,'a');
            switch upper(pv.mode)
                case "MOVIE"
                    % Show the movie with the decoded pupil on top.
                    fprintf('Opening movie file...\n')
                    movie = VideoReader(fullfile(fldr,cTpl.filename)); %#ok<TNMLP>
                    fprintf('done.\n Press Ctrl-C to stop, or close the window to move to the next movie in the table\n');
                    phi = linspace(0,2*pi,100)';
                    colormap gray
                    for frameCntr = pv.frameStart:pv.frameStep:min(pv.frameStop,cTpl.nrsamples)
                        if  ~ishandle(hFig);break;end
                        ax = gca;
                        frame = movie.read(frameCntr);
                        hold off
                        imagesc(frame);
                        hold on
                        plot(x(frameCntr),y(frameCntr),'r*');
                        radius = sqrt(a(frameCntr)/pi);
                        line(x(frameCntr)+radius.*cos(phi),y(frameCntr)+radius.*sin(phi),'Color','g')
                        xlabel 'X (pixels)';
                        ylabel 'Y (pixels)';
                        text(ax,max(ax.XLim), min(ax.YLim),sprintf('%s       %dx - Frame #%d/%d',movie.Name,pv.frameStep,frameCntr,cTpl.nrsamples),'HorizontalAlignment','Right','VerticalAlignment','top','Color','y','FontWeight','Bold','FontSize',12,'Interpreter','none')                     
                        drawnow;
                    end
                case "TRAJECTORY"
                    % Show trajectory x,y, area on the screen.
                    scatter(x,y,a,'ko');
                    set(gca,'XLim',[1 cTpl.info.width],'Ylim',[1 cTpl.info.height]);
                    xlabel 'X (pixels)';
                    ylabel 'Y (pixels)';
                    title(sprintf('mean quality= %.2f  NaN-Frac=%.2f',mean(q,'omitnan'),mean(isnan(a))));
                case "TIMECOURSE"
                    % Show x,y, area as a function of time.
                    T=tiledlayout(3,1,"TileSpacing","tight");
                    t =linspace(cTpl.time(1),cTpl.time(end),cTpl.time(3));
                    nexttile(T)
                    plot(t,x/ cTpl.info.width);
                    hold on
                    plot(t,y/ cTpl.info.height);
                    ylim([0 1]);
                    ylabel 'Position (frac)'
                    legend('x','y')
                    nexttile(T)
                    scaledA = (a-mean(a,"omitnan"))/std(a,0,"omitnan");
                    plot(t,scaledA);
                    ylabel 'Area (z-score)'
                    legend('area')
                    nexttile(T)

                    plot(t,q);
                    ylabel 'Position Quality ([0 1])'
                    xlabel 'Time (s)'
                    ylim([0 1]);
                    legend('quality')
                    title(T,sprintf('mean quality= %.2f  NaN-Frac=%.2f',mean(q,'all','omitnan'),mean(isnan(a),"all")));
            end
        end
    else
        fprintf('No C data for %s. Populate ns.C first? \n',exptName)
    end  % has prepped data in ns.C
end % Experiment
end

function v = channelByName(tpl,name)

ix = strcmpi({tpl.name},name);
v = tpl(ix).signal;
end

