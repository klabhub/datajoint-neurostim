%{
# Stores eye position information 
-> ns.Experiment # Corresponding experiment
-> sbx.EyeParms  # The parameters that define the pose extraction process.
---
x :longblob  # The x position of the pupil center; [nrTimePoints 1]
y :longblob  # The y position of the pupil center; [nrTimePoints 1]
a: longblob  # The pupil area [nrTimePoints 1]
quality : longblob # The quality of the estimation at each time point [nrTimePoints 1]
manualqc = NULL : smallint # quantify the overall quality based on manual inspection. 
nrtimepoints :  int unsigned # Number of time points in the pose estimation
width = NULL : float # Width of the camera image
height = NULL  : float # Height of the camera image
framerate = NULL : float # Framerate of the movie
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
            mvFile =  sbx.Eye.movieFile(key);
            if ~exist(mvFile,"file")
                error('%s file not found. (Is NS_ROOT set correctly (%s)?)',mvFile,getenv('NS_ROOT'));
            end
            movie = VideoReader(mvFile);
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
            switch upper(key.tag)
                case 'IMFINDCIRCLES'
                    %% Pupil tracking, using imfindcircles
                    movie = openMovie(tbl,key);
                    [x,y,a,quality,nrT,w,h,fr] = sbx.Eye.imfindcircles(movie, parms);
                otherwise
                    if startsWith(key.tag,'DLC','IgnoreCase',true)
                        % Any tag that starts with DLC is processed with
                        % DLC to allow DLC model comparisons.
                        mvFile =  sbx.Eye.movieFile(key);
                        [x,y,a,quality,nrT] = sbx.Eye.dlc(mvFile, parms);
                        w=NaN;h=NaN;fr=NaN;
                    else
                        error('Unknown %d tag',key.tag);
                    end
            end
            tpl = mergestruct(key,struct('x',x,'y',y,'a',a,'quality',quality,'nrtimepoints',nrT,'width',w,'height',h,'framerate',fr));
            insert(tbl,tpl);
        end
    end

    methods (Static)
        function v= movieFile(key)
            % It is possible that there is more than one _eye file; pick the smallest one
            % (presumably this is a preprocess/compressed version)
            fldr = folder(ns.Experiment& key);
            minSize = Inf;
            v='';
            for f=fetch(ns.File & key & 'filename LIKE ''%_eye%''','filename')'
                ff =fullfile(fldr,f.filename);
                if exist(ff,"file")
                    d = dir(ff);
                    if d.bytes<minSize
                        minSize= d.bytes;
                        v= ff;
                    end
                end
            end
        end
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


        function [x,y,a,quality,nrFrames] =dlc(mvFile,parms)
            % Use DeepLabCut to determine the pupil position. The parms
            % (from sbx.EyeParms) must specify the following parameters of
            % the analyze_videos function in DLC:
            %    .config
            %    .shuffle (1,1) double {mustBeNonnegative,mustBeInteger}
            %    .trainingsetindex (1,1) double {mustBeNonnegative,mustBeInteger}
            %    .TFGPUinference  (1,1) logical
            % All other analyze_videos parameters use the default value and
            % gputouse is determined on the fly.
            %
            % The parms struct should also specify the suffix it expects
            % the DLC output to have. This is used to determine which
            % DLC output file to use. This will looks something like this
            %
            % 'DLC_resnet50_EyeTrackerSep27shuffle1_650000'
            %
            % If the remote cluster uses singularity to run DLC , specify the full
            % path to the singularity (SIF) file in parms.singularity
            %
            % If the remote cluster uses a conda environment to run DLC,
            % specify the name of the environment in parms.conda.env. If
            % conda activate env would fail on your system (because your
            % bashrc does not initialize conda, you can add a command to
            % execute as parms.conda.init (e.g. source ~/.condainit) if
            % ~/.condainit contains the initialization code that is normally in bashrc).
            %
            % Ultimately this function will call python with the DLC command constructed from
            % the parms, which will write the output to the same folder as
            % the video file.
            %
            % Matlab then reads the csv files, does some postprocessing to determine
            % pupil center and area and stores the results in the Eye table.
            %
            % I trained an pupil tracker DLC model on top,left, right, and
            % bottom points of the pupil. The current postprocessing
            % determines the intersection between top-bottom and left-right
            % lines as the pupil center, and the surface of the trapezoid
            % with these four corners as the area. The quality measure is
            % the minimum of the likelihoods that DLC assigned to each of
            % the four points.
            %
            % More advanced models for pupil tracking could be integrated
            % here by changing the postprocessing code.
            %
            arguments
                mvFile  (1,1) string
                parms   (1,1) struct
            end


            mvFile = strrep(mvFile,'\','/');
            [videoFolder,videoFile,videoType]= fileparts(mvFile);
            csvFile = fullfile(videoFolder,[videoFile + parms.suffix + ".csv"]);

            if exist(csvFile,"File")
                % Skip running DLC
                fprintf('DLC output (%s) already exists. Adding to the table.\n',csvFile);
            else

                %% Construct the Python command
                % In python we import deeplabcut, then call the analyze_videos
                % function, and then exit()
                % analyze_vidoes has many input arguments, we're specifying
                % only the ones in parms, the rest take their default value.
                %
                % analyze_videos(config, videos, videotype='', shuffle=1, trainingsetindex=0, gputouse=None, save_as_csv=False,
                % in_random_order=True - Not used as we're passing one video file at a time
                % destfolder=None - Not used so the results will be written to the same folder as the vidoes
                %  batchsize=None, Not used
                % cropping=None, Not used
                % dynamic=(False, 0.5, 10), modelprefix='',
                % robust_nframes=False,
                % allow_growth=False,
                % use_shelve=False,
                % auto_track=True,
                % n_tracks=None,
                % calibrate=False,
                % identity_only=False,
                % use_openvino=None)


                if canUseGPU
                    % In case we really need to determine which gpu is availabel, something like this may work
                    %  [~,ps] = system('ps -u bart');
                    %   ps =strsplit(ps ,{' ','\n'});
                    %   ps(cellfun(@isempty,ps)) =[];
                    %   ps = reshape((ps),4,[])';
                    %   T= table(ps(2:end,1),ps(2:end,2),ps(2:end,3),ps(2:end,4),'VariableNames',ps(1,:));
                    % And compare that with
                    % nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory --format=csv
                    gputouse = '0'; % Manual says to give the index, but '0' seems to work even if 2 is assigned to Matlab?.
                    nvoption = '--nv'; % Singularity only
                else
                    gptouse = 'none'; %#ok<NASGU>
                    nvoption = '';
                end

                videoType=extractAfter(videoType,'.');
                % The python command is always the same
                pythonCmd = sprintf("import deeplabcut;deeplabcut.analyze_videos('%s',['%s'],videotype='%s',shuffle=%d,trainingsetindex=%d,gputouse=%s,save_as_csv=1,TFGPUinference=%d);exit();",parms.config,mvFile,videoType,parms.shuffle,parms.trainingsetindex,gputouse,parms.TFGPUinference);
                % But the call to run DLC can differ:
                if isfield(parms,'singularity')
                    % Run DLC in a singularity container
                    cmd = sprintf('singularity exec %s %s python -Wdefault -c "%s"',nvoption,parms.singularity, pythonCmd);
                elseif isfield(parms,'conda')
                    % Run in a conda environment (recommended)
                    cmd = sprintf('conda activate %s; python -Wdefault -c "%s"', parms.conda.env,pythonCmd);
                    if ~isempty(parms.conda.init)
                        % Prepend cona initialization code provided in the
                        % parms
                        cmd = [parms.conda.init ';' cmd];
                    end
                else
                    % Python install without an environment.
                    cmd = sprintf('python -Wdefault -c "%s"', pythonCmd);
                end

                %% Start DLC
                try
                    fprintf('Runing system command:\n\n %s \n\n',cmd);
                    [status] = system(cmd,'-echo');
                catch me
                    fprintf('DLC failed %s\n',me.message);
                end

                if status~=0
                    fprintf('DLC returede status %d\n',status);
                end
            end



            % Read the csv file
            if exist(csvFile,"file")
                T = readdlc(csvFile);
                % The DLC model determins the left, right, top, and
                % bottom points of the pupil
                % Determine the intersection of the line from top to
                % bottom and the line from left to right to define the
                % center (the intersection) and the area (the trapezoid
                % spanned by the four points).
                slopeTopBottom = (T.topy - T.bottomy) ./ (T.topx - T.bottomx);
                intersectTopBottom = T.bottomy - slopeTopBottom .* T.bottomx;

                slopeLeftRight = (T.righty - T.lefty)  ./ (T.rightx - T.leftx);
                intersectLeftRight = T.lefty - slopeLeftRight .* T.leftx;

                % Find intersection point
                x  = (intersectLeftRight - intersectTopBottom) ./ (slopeTopBottom - slopeLeftRight);
                y = slopeTopBottom .* x+ intersectTopBottom;

                topToBottom = sqrt((T.topy - T.bottomy).^2+ (T.topx - T.bottomx).^2);
                leftToRight = sqrt((T.righty - T.lefty).^2+ (T.rightx - T.leftx).^2);
                a = topToBottom.*leftToRight;

                quality = min([T.toplikelihood  T.rightlikelihood T.leftlikelihood T.bottomlikelihood],2,'omitnan');
                nrFrames= height(T);
            else
                dir(videoFolder);
                error('The expected DLC output file (%s) was not found. Check the suffix (%s) in sbx.EyeParms',csvFile,parms.suffix);
            end


        end
    end

end
