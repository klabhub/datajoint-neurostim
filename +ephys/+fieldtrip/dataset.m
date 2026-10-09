function data = dataset(key,pv)
%dataset - Build a FieldTrip raw data structure from ns.Epoch
%   DATA = DATASET(KEY) extracts one ns.Epoch tuple and returns an
%   epoched FieldTrip raw data structure.
%
%   DATA = DATASET(KEY,dimension=VALUE) additionally restricts the epoch
%   to the specified dimension.
%
%   The returned structure contains one matrix in DATA.trial and one time
%   vector in DATA.time for each retained trial. DATA.trialinfo contains
%   the original trial number, alignment onset, and numeric condition code.
%   DATA.epoch contains the corresponding trial, onset, and condition
%   metadata.
%
%   See also ns.Epoch, ns.EpochChannel, ft_datatype_raw

arguments
    key (1,1)
    pv.dimension (1,1) string = ""
end

epochKey = key;
if pv.dimension ~= ""
    epochKey.dimension = char(pv.dimension);
end
epochKey = ns.stripToPrimary(ns.Epoch, epochKey);

epochRel = ns.Epoch & epochKey;
nrEpochRows = count(epochRel);
assert(nrEpochRows == 1, ...
    "Expected one ns.Epoch tuple; found %d.", nrEpochRows);

epochTime = fetch1(epochRel, 'time');
assert(isnumeric(epochTime) && numel(epochTime) == 3 && ...
    epochTime(3) == round(epochTime(3)) && epochTime(3) >= 1, ...
    'ns.Epoch.time must be [start stop nrSamples].');

nrSamples = epochTime(3);
time = linspace(epochTime(1), epochTime(2), nrSamples);

cRel = ns.C & ns.stripToPrimary(ns.C, epochKey);

if nrSamples > 1 && time(2) ~= time(1)
    fsample = 1 / (time(2) - time(1));
else
    fsample = cRel.samplingRate;
end

continuousTime = double(cRel.time(:))/1000;% Convert to s
assert(~isempty(continuousTime) && all(isfinite(continuousTime)) && ...
    all(diff(continuousTime) > 0), ...
    'The continuous ns.C time axis must be finite and strictly increasing.');

[channel,trial,signal,onset] = fetchn(...
    ns.EpochChannel & epochRel, 'channel', 'trial', 'signal', 'onset');
assert(~isempty(signal), 'The selected ns.Epoch contains no data.');

channel = double(channel(:));
trial = double(trial(:));
onset = double(onset(:))/1000; % Convert to seconds for FT
if ~iscell(signal)
    signal = num2cell(signal,2);
end
signal = signal(:);

assert(numel(channel) == numel(trial) && ...
    numel(trial) == numel(signal) && numel(signal) == numel(onset), ...
    'EpochChannel fields have inconsistent row counts.');

trialValues = unique(trial, 'stable');
channelValues = unique(channel, 'stable');
nrTrials = numel(trialValues);
nrChannels = numel(channelValues);

trialIndex = zeros(size(trial));
channelIndex = zeros(size(channel));
[~,trialIndex(:)] = ismember(trial, trialValues);
[~,channelIndex(:)] = ismember(channel, channelValues);

dataTrial = cell(1,nrTrials);
seen = false(nrChannels,nrTrials);
for iTrial = 1:nrTrials
    dataTrial{iTrial} = zeros(nrChannels,nrSamples,'like',signal{1});
end

for iRow = 1:numel(signal)
    thisSignal = signal{iRow};
    assert(isvector(thisSignal) && numel(thisSignal) == nrSamples, ...
        ['EpochChannel signal for trial %d/channel %d has %d samples; ' ...
        'expected %d.'], trial(iRow), channel(iRow), ...
        numel(thisSignal), nrSamples);

    iChannel = channelIndex(iRow);
    iTrial = trialIndex(iRow);
    assert(~seen(iChannel,iTrial), ...
        'Duplicate EpochChannel row for trial %d/channel %d.', ...
        trial(iRow), channel(iRow));
    dataTrial{iTrial}(iChannel,:) = reshape(thisSignal,1,[]);
    seen(iChannel,iTrial) = true;
end
assert(all(seen,'all'), ...
    'EpochChannel does not contain a complete trial-by-channel grid.');

cChannels = fetch(ns.CChannel & cRel, 'channel', 'name', 'channelinfo');
[isPresent,cIndex] = ismember(channelValues, double([cChannels.channel]));
assert(all(isPresent), ...
    'Some epoch channels are missing from the corresponding ns.CChannel.');

labels = cell(nrChannels,1);
for iChannel = 1:nrChannels
    labels{iChannel} = channelLabel(cChannels(cIndex(iChannel)), ...
        channelValues(iChannel));
end
channelInfo = {cChannels(cIndex).channelinfo};
elec = electrodeStruct(channelInfo, labels);

trialOnset = nan(nrTrials,1);
for iTrial = 1:nrTrials
    trialOnset(iTrial) = onset(find(trial == trialValues(iTrial),1));
end

% FieldTrip sampleinfo refers to the sample numbers in the original
% continuous recording, not to the local sample numbers within each epoch.
% EpochChannel.onset is the absolute continuous-data time (in ms) at the
% alignment event, whereas epochTime is relative to that event.
sampleinfo = zeros(nrTrials,2);
for iTrial = 1:nrTrials
    absoluteWindow = trialOnset(iTrial) + epochTime(1:2);
    [~, sampleinfo(iTrial,1)] = min(abs(continuousTime - absoluteWindow(1)));
    [~, sampleinfo(iTrial,2)] = min(abs(continuousTime - absoluteWindow(2)));
    assert(sampleinfo(iTrial,1) <= sampleinfo(iTrial,2), ...
        'Epoch %d maps to an invalid continuous sample interval.', ...
        trialValues(iTrial));
end

condition = strings(nrTrials,1);
conditionRows = fetch(ns.DimensionTrial & epochRel, 'trial', 'name');
for iRow = 1:numel(conditionRows)
    iTrial = find(trialValues == double(conditionRows(iRow).trial),1);
    if ~isempty(iTrial)
        condition(iTrial) = string(conditionRows(iRow).name);
    end
end

conditionValues = unique(condition(condition ~= ""), 'stable');
conditionCode = nan(nrTrials,1);
for iCondition = 1:numel(conditionValues)
    conditionCode(condition == conditionValues(iCondition)) = iCondition;
end

epoch = repmat(struct('trial',[],'onset',[],'condition',""),nrTrials,1);
for iTrial = 1:nrTrials
    epoch(iTrial).trial = trialValues(iTrial);
    epoch(iTrial).onset = trialOnset(iTrial);
    epoch(iTrial).condition = condition(iTrial);
end

data = struct();
data.label = labels;
data.elec = elec;
data.trial = dataTrial;
data.time = repmat({time},1,nrTrials);
data.fsample = fsample;
data.trialinfo = [trialValues(:), trialOnset, conditionCode];
data.trialinfo_label = {'trial','onset','condition'};
data.epoch = epoch;
data.sampleinfo = sampleinfo;
data.hdr = struct('Fs',fsample,'nChans',nrChannels, ...
    'nSamples',nrSamples,'nTrials',nrTrials,'label',{labels});
data.cfg = struct();
data.cfg.dataset = sprintf('%s/%s/%s/%s/%s',epochKey.subject,epochKey.session_date,epochKey.starttime,epochKey.etag,condition);
data.cfg.ns = struct('key',epochKey,'epoch',epochKey, ...
    'trial',trialValues,'condition',condition);
end

function label = channelLabel(channelRow,channelNumber)
% Return the best available channel label from ns.CChannel metadata.
if isfield(channelRow,'name') && ...
        strlength(string(channelRow.name)) > 0
    label = char(string(channelRow.name));
    return
end

info = channelRow.channelinfo;
if isstruct(info)
    if isfield(info,'labels') && ~isempty(info.labels)
        value = info.labels;
        if iscell(value)
            value = value{1};
        end
        label = char(string(value));
        return
    elseif isfield(info,'label') && ~isempty(info.label)
        value = info.label;
        if iscell(value)
            value = value{1};
        end
        label = char(string(value));
        return
    end
end

label = sprintf('chan%d',channelNumber);
end

function elec = electrodeStruct(channelInfo, labels)
% Convert stored channel metadata into a FieldTrip electrode structure.
nrChannels = numel(labels);
position = nan(nrChannels,3);

for iChannel = 1:nrChannels
    info = channelInfo{iChannel};
    if isfield(info,'X') && isfield(info,'Y') && isfield(info,'Z')
        position(iChannel,:) = [double(info.X), double(info.Y), double(info.Z)];
    elseif isfield(info,'chanpos')
        position(iChannel,:) = reshape(double(info.chanpos),1,3);
    elseif isfield(info,'elecpos')
        position(iChannel,:) = reshape(double(info.elecpos),1,3);
    elseif isfield(info,'pos')
        position(iChannel,:) = reshape(double(info.pos),1,3);
    end
end

assert(all(isfinite(position),'all'), ...
    'Electrode position information is incomplete in ns.CChannel.');

elec = struct();
elec.label = labels;
elec.chanpos = position;
elec.elecpos = position;
elec.pnt = position;
elec.channelinfo = channelInfo;
end
