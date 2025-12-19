function [signal,time,channelInfo,recordingInfo] = readMovie(key,parms)
% Read from a movie files to create an entry in ns.C that stores the
% mean intensity per frame. The main goal is to determine nsTime for each
% movie frame so that it can be mapped to a trial with ns.C.align()
% 
% The first column of signal = 1:nrFrames, the second is the mean intensity
% of each frame.
%
% The parms struct could be used to specify additional preprocessing, but
% it not used currently.
%
% BK - Oct 2023
arguments
    key % The key (ns.File and ns.CParm tuple)
    parms (1,1) struct  %#ok<INUSA> % The preprocessing parameters
end

%% Fetch and open the file (ns.C has already checked that it exists)
filename = fullfile(folder(ns.Experiment &key),fetch1(ns.File &key,'filename'));
tic
fprintf('Reading movie file %s. ', filename)
movie = VideoReader(filename); %
frames= read(movie);
frameDim = ndims(frames);
nrFrames= size(frames,frameDim);
fprintf('Done in %d seconds.\n ',round(toc))
%% 
% Create a signal that represenst the frame counter and the mean intensity
signal = [(1:nrFrames)' squeeze(mean(frames,1:(frameDim-1),"omitnan"))];
channelInfo =struct('name',{'frame','intensity'}','nr',{1,2}');

%% Map frames to nsTime 
if exists(ns.Plugin & key & 'plugin_name=''camera''')
    % A movie captured by the neurostim camera plugin
    prms = get(ns.Experiment & key,'camera');
    withData= cellfun(@numel,prms.camera.firstvideoframe);
    %% This code is for SaveDuringTrial mode : 
    % one trigger before the trial starts triggers collection of a fixed number of frames.
    % But the afterFrame callbacks do not always run (hence nans and
    % negative offsets) even when the movie contains a full set of frames.
    % The timing of frames is probably still fine (as long as no subframe
    % accuracy is expected).
    offsets = cat(1,prms.camera.firstvideoframe{withData>0})';
    [nrFramesPerTrial,~] =size(offsets);
    out= offsets<=0 | isnan(offsets);

    mFrameDuration = mean(offsets(~out));
    sdFrameDuration  = std(offsets(~out));
    z = sdFrameDuration/mFrameDuration;
    fprintf('%.0f%% of frames have NaN offsets, z-score of the variation in frameduration is %.2f\n',100*mean(isnan(offsets),"all"),z)
    actualFrameRate = 1000./mFrameDuration;
    % Assuming a constant framerate 
    time = prms.camera.firstvideoframeNsTime(withData>0)'+ (0:nrFramesPerTrial-1)'*mFrameDuration;
else
    % Assume movie started with the first frame in the first trial
    error('No camera plugin found. Alignment to Neurostim clock not possible?');
end


recordingInfo = struct('nrFrames',movie.NumFrames,'videoFormat',movie.VideoFormat,'width',movie.Width,'height',movie.height,'framerate',actualFrameRate);
% Regular sampling:  reduce time representation
time = [time(1) time(end) nrFrames];
% Reduce storage (ns.C.align converts back to double)
signal  = single(signal);

end



