function [signal,time,channelInfo,recordingInfo] = readMovie(key,parms)
% Read from a SBX movie files (_eye and _ball) to create an entry in ns.C
% that stores the eye position and ball speed.
%
% The parms structures contain the preprocessing parameters 
% 
%% BALL:
%   .method = 'xcorr' or 'phasecorr'
%   .scaleFactor = By how much to scale down (imrresize) the movie.
%   .minPixels = Scaling always leaves this number of pixels in width and
%                   height of the movie.
% 
%% EYE:
%   .method = 'imfindcircles', 'dlcXXX'
% For method = imfindcircles, specify the options of the Matlab imfindcircles
% function (listed are the defaults):
%    .objectPolarity  ['bright']
%    .method          ['PhaseCode']
%    .sensitivity    [0.85]
%    .edgeThreshold  []
%    .radiusRange   [6 20]
%
% To use a DeepLabCut trained model, use a method that starts with dlc.
% The parms struct determines how DLC is run (e.g. using conda) and which
% model to use. For details, see sbxDlc in sbx.readMovie
%
% See also sbxXcorr,sbxPhasecorr,sbxImfindcircles, sbxDlc in sbx.readMovie
%
% BK - Oct 2023
arguments
    key % The key (ns.File and ns.CParm tuple)
    parms (1,1) struct  % The preprocessing parameters
end

%% Fetch and open the file (ns.C has already checked that it exists)
fldr = folder(ns.Experiment &key);
filename = fullfile(fldr,fetch1(ns.File &key,'filename'));
tic
movie = VideoReader(filename); %
fprintf('Done in %d seconds.\n ',round(toc))


%% Map frames to nsTime
% SBX movies are matched to the frames of the TPI acquisition. Those are
% registered on the neurostim clock by the mdaq plugin. Check if they are
% already stored as events in the ns.PluginParameter, if not, read the bin
% file to exract.
if exists(ns.PluginParameter & key & 'plugin_name=''mdaq''' & 'property_name=''laserOn''')
    % Digital events have been created to represent laserOn triggers
    time = get(ns.Experiment &key,'mdaq','prm','laserOn','what','clockTime');
else
    % Read the bin file
    nsFile  = file(ns.Experiment& key);
    binFile =strrep(nsFile,'.mat','.bin');
    if exist(binFile,"file")
        fprintf('Reading bin data file %s. ', binFile)
        tic
        load(nsFile,'c');
        T=  readBin(c.mdaq,file=binFile);
        fprintf('Done in %d seconds.\n ',round(toc))
        % Determine when the laser triggers occcured
        time = T.nsTime([false; diff(T.laserOnDig)>0]);
        if contains(filename,'_ball')
            % (Ball triggers on both up and down)
            time = [time; T.nsTime([false; diff(T.laserOnDig)<0])];
            time = sort(time);
        end
        nrTriggers= numel(time);
    else
        error('No mdaq .bin file found to process %s movie',filename);
    end
end

% Sanity checks
SLACKFRAMES = 2; % Allow this many missing frames
assert(abs(nrTriggers - movie.NumFrames)<SLACKFRAMES,"Movie %s has %d frames, but %d triggers ",filename,movie.NumFrames,nrTriggers);
dt = seconds(diff(time));
mFrameDuration = mean(dt);
sdFrameDuration  = std(dt);
z = sdFrameDuration/mFrameDuration;
fprintf('Z-score of the variation in frameduration is %.2f\n',z)
actualFrameRate = 1000./mFrameDuration;
recordingInfo = struct('nrFrames',movie.NumFrames,'videoFormat',movie.VideoFormat,'width',movie.Width,'height',movie.height,'framerate',actualFrameRate);

%% Analyze the movies using functions defined below
if contains(filename,'_ball')
    switch upper(parms.method)
        case 'XCORR'
            % cross correaltion
            [velocity,quality] =sbxXcorr(movie,parms);
        case 'PHASECORR'
            % phase correlation
            [velocity,quality] =sbxPhasecorr(movie,parms);
        otherwise            
            error('Unknown method %s ',parms.method);
    end
    % Prep for ns.C
    signal = [velocity quality];
    channelInfo= struct('name',{'velocity','quality'},'nr',{1,2});
elseif contains(filename,'_eye')
    switch upper(parms.method)
        case 'IMFINDCIRCLES'
            %% Pupil tracking, using imfindcircles
            [x,y,a,quality] = sbxImfindcircles(movie, parms);
        otherwise
            if startsWith(parms.method,'DLC','IgnoreCase',true)
                % Any method that starts with DLC is processed with
                % DLC to allow DLC model comparisons.
                [x,y,a,quality] =sbxDlc(mvFile, parms);
            else
                error('Unknown method %s ',parms.method);
            end
    end
    signal =  [x y a quality];
    channelInfo= struct('name',{'x','y','a','quality'},'nr',{1,2,3,4});
else
    error('Unknown SBX movie file %s',filename)
end



%% Package to add to ns.C
% Regular sampling:  reduce time representation
time = [time(1) time(end) nrFrames];
% Reduce storage (ns.C.align converts back to double)
signal  = single(signal);
end


%% Ball Movie Analysis functions
function [velocity,quality] =sbxPhasecorr(movie,parms)
% Use a phase correlation algorithm to determine the Ball velocity
arguments
    movie (1,1) VideoReader
    parms (1,1) struct
end
useGPU = canUseGPU;  % If a GPU is available we'll use it.

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

fprintf('Ball tracking phasecorr analysis (useGPU: %d)\n',useGPU);
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

function [velocity,quality] =sbxXcorr(movie,parms)
% Use cross-correlation to determine velocity,
arguments
    movie (1,1) VideoReader
    parms (1,1) struct
end
useGPU = canUseGPU;  % If a GPU is available we'll use it.

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

fprintf('Ball tracking analysis (useGPU: %d)\n',useGPU);
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


%% Eye Movie Analysis
function [x,y,a,quality] =sbxImfindcircles(movie,parms)
% The imfindcircles tool defines the following parameters,
% which should be specified in the parms struct.  These are the
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
    parms (1,1) struct
end
%% Use parpool if available
if isempty(gcp('nocreate'))
    nrWorkers = 0;
else
    nrWorkers  = gcp('nocreate').NumWorkers;
end

nrFrames = movie.NumFrames;
%% Initialize output vars
x       = nan(nrFrames,1);
y       = nan(nrFrames,1);
a       = nan(nrFrames,1);
quality  = nan(nrFrames,1);

fprintf('Eye tracking analysis on %d workers\n',nrWorkers)
% Read the movie into memory (could check that we have
% enough...)
frames= single(movie.read([1 nrFrames]));
parfor (f=1:nrFrames,nrWorkers)
    %for f=1:nrT  % debug
    [center,radius,thisQuality] = sbxImfindcircles(frames(:,:,f),parms.radiusRange,'ObjectPolarity',parms.objectPolarity,'Method',parms.method,'Sensitivity',parms.sensitivity,'EdgeThreshold',parms.edgeThreshold); %#ok<PFBNS>
    if ~isempty(center)
        [quality(f),idx] = max(thisQuality); % pick the circle with best score
        x(f) = center(idx,1);
        y(f) = center(idx,2);
        a(f) = pi*radius(idx)^2;
    end
end
end


function [x,y,a,quality] =sbxDlc(movie,parms)
% Use DeepLabCut to determine the pupil position. The parms
% must specify the following parameters of the analyze_videos function in DLC:
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
    movie  (1,1) VideoReader
    parms   (1,1) struct
end

mvFile = movie.File;
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
    [x,y,a,quality] = postprocessPupilTracker(T,bodyparts);
else
    dir(videoFolder);
    error('The expected DLC output file (%s) was not found. Check the suffix (%s) in the parms. ',csvFile,parms.suffix);
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
