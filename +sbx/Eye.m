%{
# Stores eye position information 
-> ns.Movie # Corresponding movie
-> sbx.EyeParms  # The parameters that define the pose extraction process.
---
x :longblob  # The x position of the pupil center; [nrTimePoints 1]
y :longblob  # The y position of the pupil center; [nrTimePoints 1]
a: longblob  # The pupil area [nrTimePoints 1]
quality : longblob # The quality of the estimation at each time point [nrTimePoints 1]
%}
%
% BK - Sept 2023.

classdef Eye < dj.Computed
    properties (Dependent)
        keySource
    end


    methods
        function v= get.keySource(~)
            % Even though the key source is .mj2, if there is an .avi that
            % is smaller, that will be used to do the processing  due to
            % the call of ns.Movie/file with -1. (smallest).
            v = (ns.Movie & 'filename LIKE ''%_eye.mj2''')*sbx.EyeParms;
        end
    end

    methods (Access=public)
        function plot(tbl,pv)
            arguments
                tbl (1,1) sbx.Eye
                pv.mode (1,1) string {mustBeMember(pv.mode,["MOVIE","TRAJECTORY","TIMECOURSE"])} = "TRAJECTORY"
            end

            for tpl = tbl.fetch('*')'
                movieParms = fetch(ns.Movie&tpl,'*');
                figName= sprintf('Eye: #%s on %s@%s',tpl.subject, tpl.session_date,tpl.starttime);
                hFig = figByName(figName);
                colormap gray
                clf

                switch upper(pv.mode)
                    case "MOVIE"
                        % Show the movie with the decoded pupil on top.
                        movie = open(ns.Movie& tpl,"smallest");
                        frameCntr = 0;
                        phi = linspace(0,2*pi,100);
                        while (movie.hasFrame) && ishandle(hFig)
                            frameCntr = frameCntr+1;
                            frame = movie.readFrame;
                            hold off
                            imagesc(frame);
                            hold on

                            x = tpl.x(frameCntr,:)';
                            y = tpl.y(frameCntr,:)';

                            plot(x,y,'r*');


                           
                                radius = sqrt(tpl.a(frameCntr)/pi);
                                line(tpl.x(frameCntr)+radius.*cos(phi),tpl.y(frameCntr)+radius.*sin(phi),'Color','g')
                           
                            xlabel 'X (pixels)';
                            ylabel 'Y (pixels)';

                            drawnow;
                        end
                    case "TRAJECTORY"
                        % Show trajectory x,y, area on the screen.
                        scatter(tpl.x, tpl.y,tpl.a,'ko');
                        set(gca,'XLim',[1 movieParms.width],'Ylim',[1 movieParms.height]);
                        xlabel 'X (pixels)';
                        ylabel 'Y (pixels)';

                        title(sprintf('mean quality= %.2f  NaN-Frac=%.2f',mean(tpl.quality,'omitnan'),mean(isnan(tpl.a))));

                    case "TIMECOURSE"
                        % Show x,y, area as a function of time.
                        T=tiledlayout(3,1,"TileSpacing","tight");
                        t =(0:movieParms.nrframes-1)/movieParms.framerate;
                        nexttile(T)
                        plot(t,tpl.x./movieParms.width);
                        hold on
                        plot(t,tpl.y/movieParms.height);
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
                        title(T,sprintf('mean quality= %.2f  NaN-Frac=%.2f',mean(tpl.quality,'all','omitnan'),mean(isnan(tpl.a),"all")));
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
                    movie = open(ns.Movie & key,"smallest");
                    [x,y,a,quality] = sbx.Eye.imfindcircles(movie, parms);
                otherwise
                    if startsWith(key.tag,'DLC','IgnoreCase',true)
                        % Any tag that starts with DLC is processed with
                        % DLC to allow DLC model comparisons.
                        mvFile =  file(ns.Movie & key,"smallest");
                        [x,y,a,quality] = sbx.Eye.dlc(mvFile, parms);
                    else
                        error('Unknown %d tag',key.tag);
                    end
            end
            tpl = mergestruct(key,struct('x',x,'y',y,'a',a,'quality',quality));
            insert(tbl,tpl);
        end
    end

    methods (Static)

        function [x,y,a,quality] =imfindcircles(movie,pv)
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


        function [x,y,a,quality] =dlc(mvFile,parms)
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
            %
            % If the remote cluster uses a conda environment to run DLC,
            % specify the name of the environment in parms.conda.env. If
            % conda activate env would fail on your system (because your
            % bashrc does not initialize conda, you can add a command to
            % execute as parms.conda.init (e.g. source ~/.condainit) if
            % ~/.condainit contains the initialization code that is normally in bashrc).
            % If the cluster does not use python, set parms.conda =""
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
            videoType=extractAfter(videoType,'.');
            csvFile = fullfile(videoFolder,videoFile + parms.suffix + ".csv");

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
            else
                gptouse = 'none'; %#ok<NASGU>
            end



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
                pythonCmd = sprintf("import deeplabcut;deeplabcut.analyze_videos('%s',['%s'],videotype='%s',shuffle=%d,trainingsetindex=%d,gputouse=%s,save_as_csv=1,TFGPUinference=%d);exit();",parms.config,mvFile,videoType,parms.shuffle,parms.trainingsetindex,gputouse,parms.TFGPUinference);
                rundlc(pythonCmd,"condaEnv",parms.conda.env,"condaInit",parms.conda.init);
            end

            if isfield(parms,'filter')
                % Generate the filtered predictions.
                % parms.filter can have the following fields
                % parms.filter.filtertype,
                % parms.filter.windowlength,
                % parms.filter.p_bound,
                % parms.filter.ARdegree,
                % parms.filter.MAdegree,
                % parms.filter.alpha
                % Missing fields will get the default values in DLC
                % deeplabcut.post_processing.filtering.filterpredictions(config, video, videotype='', shuffle=1, trainingsetindex=0, filtertype='median', windowlength=5, p_bound=0.001, ARdegree=3, MAdegree=1, alpha=0.01, save_as_csv=True, destfolder=None, modelprefix='', track_method='')

                csvFile = fullfile(videoFolder,videoFile + parms.suffix + "_filtered.csv");
                if exist(csvFile,"file")
                    fprintf('Filtered DLC output (%s) already exists. Adding to the table.\n',csvFile);
                else
                    pythonCmd = sprintf("import deeplabcut;deeplabcut.filterpredictions('%s',['%s'],videotype='%s',shuffle=%d,trainingsetindex=%d,save_as_csv=1" ,parms.config,mvFile,videoType,parms.shuffle,parms.trainingsetindex);
                    fn = fieldnames(parms.filter);
                    filterArgs = cell(1,numel(fn));
                    for i=1:numel(fn)
                        value =parms.filter.(fn{i});
                        if isnumeric(value)
                            value=num2str(value);
                        else
                            value = ['''' value '''']; %#ok<AGROW>
                        end
                        filterArgs{i}  = sprintf('%s=%s',fn{i},value);
                    end
                    pythonCmd = pythonCmd +"," + strjoin(filterArgs,",")+ ");";
                    rundlc(pythonCmd,"condaEnv",parms.conda.env,"condaInit",parms.conda.init);
                end
            end

            % Read the csv file
            if exist(csvFile,"file")
                [T,bodyparts] = readdlc(csvFile);
                [x,y,a,quality] = sbx.Eye.postprocessPupilTracker(T,bodyparts);
            else
                dir(videoFolder);
                error('The expected DLC output file (%s) was not found. Check the suffix (%s) in sbx.EyeParms',csvFile,parms.suffix);
            end
        end



        function [x,y,a,quality] = postprocessPupilTracker(T,bodyparts)
            % Given a table from a csv file output created by a specific
            % DLC model that determins the left,
            % right, top, and bottom points of the pupil, determine the intersection of the line from top to
            % bottom and the line from left to right to define the
            % center (the intersection) and the area (the trapezoid
            % spanned by the four points).
            isTop = strcmpi(bodyparts,'top');
            isBottom = strcmpi(bodyparts,'bottom');
            isLeft = strcmpi(bodyparts,'left');
            isRight = strcmpi(bodyparts,'right');

            slopeTopBottom = (T.y(:,isTop) - T.y(:,isBottom)) ./ (T.x(:,isTop) - T.x(:,isBottom));
            intersectTopBottom = T.y(:,isBottom) - slopeTopBottom .* T.x(:,isBottom);

            slopeLeftRight = (T.y(:,isRight)- T.y(:,isLeft))  ./ (T.x(:,isRight) - T.x(:,isLeft));
            intersectLeftRight = T.y(:,isLeft) - slopeLeftRight .* T.x(:,isLeft);

            % Find intersection point
            x  = (intersectLeftRight - intersectTopBottom) ./ (slopeTopBottom - slopeLeftRight);
            y = slopeTopBottom .* x+ intersectTopBottom;
            topToBottom = sqrt((T.y(:,isTop) - T.y(:,isBottom)).^2+ (T.x(:,isTop) - T.x(:,isBottom)).^2);
            leftToRight = sqrt((T.y(:,isRight) - T.y(:,isLeft)).^2+ (T.x(:,isRight) - T.x(:,isLeft)).^2);
            a = topToBottom.*leftToRight;
            quality = min(T.likelihood,[],2,'omitnan');           
        end

    end

end