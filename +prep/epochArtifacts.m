function [signal,badBy,details] = epochArtifacts(signal,srate,options)
%epochArtifacts - Reject epochs or mask artifact windows
%   [SIGNAL,BADBY,DETAILS] = epochArtifacts(SIGNAL,SRATE,OPTIONS) accepts
%   samples-by-trials-by-channels data and the artifact options stored in
%   ns.EpochParm. The default action "reject" returns trial flags in BADBY;
%   the caller removes those trials. This function retains array dimensions.
%
%   OPTIONS.action="nan" assesses temporary windows together using
%   prep.artifactDetection and masks flagged intervals across all channels.
%   OPTIONS.window, step, and padding are milliseconds, defaulting to 1000,
%   500, and 100. The step must not exceed the window. The final window is
%   included; an epoch shorter than the window uses its available samples.
%
%   OPTIONS.epoch_no maps array columns to original trial numbers;
%   OPTIONS.exclude skips those original trials. Windows containing missing
%   or infinite samples are masked without entering the detector statistics.
%   DETAILS contains window flags, merged sample intervals, and masked
%   fractions in trial order. BADBY remains empty for action "nan".
%
%   See also prep.artifactDetection, ns.Epoch
arguments
    signal (:,:,:) {mustBeFloat,mustBeReal}
    srate (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBePositive}
    options (1,1) struct
end

defaults = struct('action',"reject",'window',1000,'step',500,'padding',100);
for name = string(fieldnames(defaults))'
    if ~isfield(options,name)
        options.(name) = defaults.(name);
    end
end
mustBeMember(string(options.action),["reject" "nan"]);
mustBeTextScalar(options.action);
validateattributes(options.window,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(options.step,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(options.padding,{'numeric'},{'scalar','real','finite','nonnegative'});
assert(options.step <= options.window,'prep:epochArtifacts:Step', ...
    'Artifact step must not exceed the detection window.');
detector = rmfield(options,{'action','window','step','padding'});
if isfield(detector,'ica')
    detector = rmfield(detector,'ica');
end
badBy = prep.badBy(detector);
details = struct('action',string(options.action),'parameters',options);
if string(options.action) == "reject"
    args = namedargs2cell(detector);
    badBy = prep.artifactDetection(permute(signal,[2 3 1]),srate,args{:});
    return
end

[nrSamples,nrTrials,nrChannels] = size(signal);
trialNumbers = (1:nrTrials)';
if isfield(detector,'epoch_no') && ~isempty(detector.epoch_no)
    trialNumbers = detector.epoch_no(:);
end
assert(numel(trialNumbers) == nrTrials && ...
    numel(unique(trialNumbers)) == nrTrials, ...
    'prep:epochArtifacts:TrialNumbers','epoch_no must identify each trial uniquely.');
excluded = false(nrTrials,1);
if isfield(detector,'exclude')
    excluded = ismember(trialNumbers,detector.exclude);
end
mask = false(nrSamples,nrTrials);
windows = table(zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'trial','trialIndex','firstSample','lastSample'});
windowFlags = prep.badBy(detector);
enabled = ~isfield(detector,'enable') || detector.enable;
if enabled && nrTrials > 0
    width = min(nrSamples,max(2,round(options.window*srate/1000)));
    assert(width >= 2,'prep:epochArtifacts:SignalLength', ...
        'Windowed detection requires at least two samples.');
    step = max(1,round(options.step*srate/1000));
    starts = unique([1:step:(nrSamples-width+1),nrSamples-width+1]);
    selectedTrials = find(~excluded);
    nrWindows = numel(starts)*numel(selectedTrials);
    windowData = zeros(nrWindows,nrChannels,width,'like',signal);
    windowTrial = repelem(selectedTrials,numel(starts));
    windowTrial = windowTrial(:);
    first = repmat(starts(:),numel(selectedTrials),1);
    last = first+width-1;
    windows = table(trialNumbers(windowTrial),windowTrial,first,last, ...
        'VariableNames',windows.Properties.VariableNames);
    for k = 1:nrWindows
        windowData(k,:,:) = permute(signal(first(k):last(k),windowTrial(k),:),[2 3 1]);
    end
    valid = all(isfinite(windowData),[2 3]);
    valid = valid(:);
    detector.epoch_no = find(valid);
    detector.exclude = [];
    if any(valid)
        args = namedargs2cell(detector);
        windowFlags = prep.artifactDetection(windowData(valid,:,:),srate,args{:});
    end
    windowFlags.missing = find(~valid);
    padding = ceil(options.padding*srate/1000);
    for k = windowFlags.all(:)'
        mask(max(1,first(k)-padding):min(nrSamples,last(k)+padding),windowTrial(k)) = true;
    end
end

intervals = zeros(0,3);
for tr = 1:nrTrials
    edges = diff([false;mask(:,tr);false]);
    first = find(edges == 1);
    last = find(edges == -1)-1;
    intervals = [intervals;repmat(trialNumbers(tr),numel(first),1),first,last]; %#ok<AGROW>
    signal(mask(:,tr),tr,:) = NaN;
end
details.windows = windows;
details.windowFlags = windowFlags;
details.intervals = array2table(intervals, ...
    'VariableNames',{'trial','firstSample','lastSample'});
details.trials = trialNumbers;
details.maskedFraction = mean(mask,1)';
end
