%{
# Stores pose information extracted by DeepLabCut
-> ns.Movie # Corresponding experiment
-> ns.PoseParm  # The parameters that define the pose extraction process.
---
pose :longblob  # The pose; [nrTimePoints nrTrackedItems]
likelihood : longblob # The quality of the estimation at each time point [nrTimePoints nrTrackedItems]
names :  blob   # The names of the tracked items (i.e. the names of the columns in pose)
%}
%
% BK - Sept 2023.

classdef Pose < dj.Computed
  
    methods (Access=public)
       
        function plot(tbl,pv)
            arguments
                tbl (1,1) sbx.Pose
                pv.mode (1,1) string {mustBeMember(pv.mode,["MOVIE","TRAJECTORY","TIMECOURSE"])} = "TRAJECTORY"
            end

            for tpl = tbl.fetch('*')'
                figName= sprintf('#%s on %s@%s',tpl.subject, tpl.session_date,tpl.starttime);
                figByName(figName);
                clf;
                if isnan(tpl.framerate)
                    tpl.framerate = 30;
                end
                if isnan(tpl.width)
                    tpl.width = max(tpl.x);
                    tpl.height =max(tpl.y);
                end

                switch upper(pv.mode)
                    case "MOVIE"
                        % Show the movie with the decoded pupil on top.
                        movie = openMovie(ns.Movie&tpl);
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

                        title(sprintf('mean quality= %.2f  NaN-Frac=%.2f',mean(tpl.quality,'omitnan'),mean(isnan(tpl.a))));

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
                        title(T,sprintf('mean quality= %.2f  NaN-Frac=%.2f',mean(tpl.quality,'all','omitnan'),mean(isnan(tpl.a),"all")));
                end
            end
        end
    end


    methods (Access = protected)
        function makeTuples(tbl,key)
            parms= fetch1(ns.PoseParms &key,'parms');
            mvFile =  movieFile(nsMovie & key);
            mvFile = strrep(mvFile,'\','/');
            [videoFolder,videoFile,videoType]= fileparts(mvFile);
            videoType=extractAfter(videoType,'.');
            expectedCsvFile = fullfile(videoFolder,videoFile + parms.suffix + ".csv");

           
            % Use DeepLabCut to determine the pose. The parms
            % (from ns.PoseParms) must specify the following parameters of
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
            % parms.suffix = 'DLC_resnet50_EyeTrackerSep27shuffle1_650000'
            %
            % parms.postprocess specifies the name of an optional function (on the
            % path) that takes a table with the DLC output (from the .csv
            % file) and does some postprocessing on the raw pose
            % parameters.
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
            

            
            if exist(expectedCsvFile,"File")
                % Skip running DLC
                fprintf('DLC output (%s) already exists. Adding to the table.\n',expectedCsvFile);
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
                else
                    gptouse = 'none'; %#ok<NASGU>
            end

                pythonCmd = sprintf("import deeplabcut;deeplabcut.analyze_videos('%s',['%s'],videotype='%s',shuffle=%d,trainingsetindex=%d,gputouse=%s,save_as_csv=1,TFGPUinference=%d);exit();",parms.config,mvFile,videoType,parms.shuffle,parms.trainingsetindex,gputouse,parms.TFGPUinference);
                rundlc(pythonCmd,"condaEnv",parms.conda.env,"condaInit",parms.conda.init);
            end

            if isfield(parms.filter)
                % Also generate the filtered predictions.
                % parms.filter can have the following fields
                % parms.filter.filtertype,
                % parms.filter.windowlength,
                % parms.filter.p_bound,
                % parms.filter.ARdegree,
                % parms.filter.MAdegree,
                % parms.filter.alpha
                % Missing fields will get the default values in DLC
                % deeplabcut.post_processing.filtering.filterpredictions(config, video, videotype='', shuffle=1, trainingsetindex=0, filtertype='median', windowlength=5, p_bound=0.001, ARdegree=3, MAdegree=1, alpha=0.01, save_as_csv=True, destfolder=None, modelprefix='', track_method='')
                pythonCmd = sprintf("import deeplabcut;deeplabcut.filterpredictions('%s',['%s'],videotype='%s',shuffle=%d,trainingsetindex=%d,save_as_csv=1" ,parms.config,mvFile,videoType,parms.shuffle,parms.trainingsetindex);
                fn = fieldnames(parms.filter);
                filterArgs = cell(1,numel(fn));
                for i=1:numel(fn)
                    value =parms.filter.(fn{i});
                    if isnumeric(value) ;value=num2str(value);end
                    filterArgs{i}  = sprintf('%s=%s',fn{i},value);
                end
                pythonCmd = pythonCmd +"," + strjoin(filterArgs,",")+ ");";
                rundlc(pythonCmd,"condaEnv",parms.conda.env,"condaInit",parms.conda.init);
            end

            % Read the csv file
            if exist(expectedCsvFile,"file")
                [T,bodyparts] = readdlc(expectedCsvFile);
                if isfield(parms,'postprocess')
                    if ~exist(parms.postprocess,"file")
                        error('Postprocessing function %s not found ',parms.postprocess);
                    end
                    % Apply user specified postprocessed
                    T = feval(parms.postprocess,T,bodyparts);
                end                
            else
                dir(videoFolder);
                error('The expected DLC output file (%s) was not found. Check the suffix (%s) in sbx.EyeParms',expectedCsvFile,parms.suffix);
            end

      
        end

        
    end

end
