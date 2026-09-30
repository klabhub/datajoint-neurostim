function EEG = dataset(key,pv)
% Uses EEGLAB plugins (mffmatlabio and fieldtrip) to read egi MFF
% files and returns an EEG struct. This is used to add EGI data to the
% database (ephys.egi.read),
% 
% Neurostim adds the following fields
% EEG.etc.neurostim.clockParms   - polyval(clockParms,EGI.time) -> converts EGI time to neurostim time
% EEG.etc.neurostim.pluginparameter -  the events stored in the MFF file  packaged as a struct array for easy insertion into the pluginparameter table (used by ephys.egi.read_
% EEG.etc.neurostim.expt = the Experiment tuple
%
% See also ephys.eeglab.addEvents, ephys.egi.read
arguments
    key (1,1)      % Experiment key or a keySource for the ns.C table
    pv.data (1,1) string % RAW, EMPTY, or a ctag
    pv.plg (1,:)  string  = string.empty % Plugin name 
    pv.prm (1,:) string   = string.empty % Event names
    pv.itag (1,1) string = "" % The ICA to load, identified by its itag. "" means no ICA will be loaded.
    pv.etag (1,1) string = "" % Epoch parameter tag. "" means return continuous data.
    pv.dimension (1,1) string = "" % Optional dimension restriction for epoch data.
    pv.readNsMeta (1,1) logical = false; % For reading the Neurostim meta data from the mff file ; required for Continuous data, optional for Epoched data.
end 

if ismember(upper(pv.data),["RAW" "EMPTY"])    
    % Ignore ctag to determine ns.C key
    key = fetch(ns.File & key & 'extension=".mff" AND NOT filename LIKE "%zcheck%.mff"','filename');
else
    key.ctag = pv.data;
    key = fetch(ns.C & key & 'filename LIKE "%.mff"','filename');    
end
assert(~isempty(key),"This experiment does not have an associated MFF file");
mffFile = fullfile(folder(ns.Experiment &key),key.filename);
mffFile= strrep(mffFile,'\','/'); % Avoid fprintf errors
assert(exist(mffFile),"MFF file %s does not exist.",mffFile); %#ok<EXIST>

% Epoch data are already segmented and preprocessed in ns.EpochChannel.
% Branch before continuous C data are fetched so epoch mode does not first
% load the complete recording into memory.
if pv.etag ~= ""
    assert(~ismember(upper(pv.data),["RAW" "EMPTY"]), ...
        'Epoch mode requires pv.data to identify an ns.C ctag.');
    EEG = epochDataset(key,pv,mffFile);
    if pv.readNsMeta
        EEG = addNeurostimMetadata(EEG,key,mffFile,[],true);
    end
    if ~isempty(pv.plg)
        EEG = ephys.egi.eeglabAddEvents(EEG,pv.plg,pv.prm);
    end
    EEG = ensureEpochAlignmentEvents(EEG);
    EEG = eeg_checkset(EEG);
    return
end


switch upper(pv.data)
    case "RAW"
        % Call mff_import directly to read everything
        eegLabSave = 0 ; % Don't save in eeglab
        correctEvents = 0; %  Don't correct events with UTF chars/
        fprintf("Using eeglab to read header and data from " +  mffFile + "...\n")
        EEG = pop_mffimport(char(mffFile),{},eegLabSave,correctEvents);
        urSrate = EEG.srate;
    otherwise        
        % Adapted code from mff_import to avoid reading the signal
        % Some pieces (that we don't currently use) are missing
        fprintf("Using eeglab to read header from " +  mffFile + "...\n")

        %  Initialize an empty standard EEGLAB structure
        EEG = eeg_emptyset();
        %  Use the ft read header v1 tool to parse metadata without reading the binary signal
        % (mff_import always reads the signal)
        mffHeader = ft_read_header(mffFile,'headerformat','egi_mff_v1');
        %  Map the XML metadata metadata onto the empty EEGLAB structure
        EEG.setname    = sprintf('%s@%sT%s',key.subject,key.session_date,key.starttime);
        EEG.srate      = mffHeader.Fs;
        EEG.nbchan     = mffHeader.nChans;
        EEG.pnts       = mffHeader.nSamples;
        EEG.trials     = 1;

        % Import channel locations from the MFF coordinates XML
        [EEG.chanlocs, EEG.ref] = mff_importcoordinates(mffFile);
        if iscell(EEG.ref)
            EEG.ref = sprintf('%s ', EEG.ref{:});
        end
        EEG.urchanlocs = EEG.chanlocs;
        if ~isempty(EEG.chanlocs)
            EEG = eeg_checkchanlocs(EEG); % put fiducials in chanfinfo
        end
        EEG=pop_chanedit(EEG, 'forcelocs',[],'nosedir','+Y');
        EEG.chaninfo.filename = 'egimff';

        % Determine when recording started
        [~, begTime]= mff_importinfo(mffFile);

        % Import the event tracks
        correctEvents=0;
        EEG.event      = mff_importevents(mffFile,begTime,EEG.srate,correctEvents);
        urSrate = EEG.srate;
        if upper(pv.data)=="EMPTY" 
           % Set the data to all-zero sparse to avoid erors on channel/time
            % selection
            EEG.data =sparse(EEG.nbchan,EEG.pnts);
        else
            %Get the data from a ctag
            cRel = ns.C & key & struct('ctag',pv.data);
            assert(exists(cRel),'No ns.C data with ctag %s found for this experiment.',pv.data);
            fprintf("Retrieving preprocessed data from ns.C (ctag=%s)\n",pv.data)
            C = fetch(ns.CChannel &cRel,'channelinfo','signal');
            EEG.chanlocs = [C.channelinfo];
            EEG.data = [C.signal]';
            [EEG.nbchan,EEG.pnts] = size(EEG.data);
            EEG.srate  = round(cRel.samplingRate);
            EEG.xmin  =0;
            EEG.xmax = (EEG.pnts-1)/EEG.srate+EEG.xmin;

            % Check if there are ICA results
            if pv.itag ~=""
                    icaKey = key;
                    icaKey.itag = pv.itag;
                    w = ns.Ica.getWeights(icaKey);  % handles both session and per-exp ICA
                    if isempty(w)
                                       fprintf("ICA with itag %s not found.\n",pv.itag);
                    elseif ~isempty(w.chanlabels)
                        % Session ICA: some channels may have been excluded
                        % during the session ICA (bad in another experiment).
                        % Find which channels are included vs excluded here.
                        allLabels = {EEG.chanlocs.labels};
                        [~, icachansind] = ismember(w.chanlabels, allLabels);
                        icachansind = icachansind(icachansind > 0);
                        assert(~isempty(icachansind), ...
                            'None of the session ICA channels found in this experiment''s chanlocs.');

                        EEG.icachansind = icachansind;
                        EEG.icasphere   = w.sphere;
                        EEG.icaweights  = w.weights;
                        EEG.icawinv     = w.winverse;
                        EEG.icaact = icaact(EEG.data(icachansind,:), ...
                        EEG.icaweights * EEG.icasphere, mean(EEG.data(icachansind,:), 2));

                        % Project excluded channels into ICA space via OLS so
                        % that pop_subcomp can clean them too.
                        % W_excl = D_excl * A' * inv(A * A')
                        % where A = icaact  [nComps x nSamps]
                        exclIdx = setdiff(1:EEG.nbchan, icachansind);
                        if ~isempty(exclIdx)
                            A  = EEG.icaact;                        % [nComps x nSamps]
                            D  = EEG.data(exclIdx, :);              % [nExcl  x nSamps]
                            % OLS via pseudo-inverse: handles rank-deficient icaact
                            % (e.g. zeroed-out components after pop_subcomp).
                            W_excl = (D * A') * pinv(A * A');       % [nExcl x nComps]
                            % Augment icawinv with projected rows for excluded channels
                            winv_aug = zeros(EEG.nbchan, size(EEG.icawinv, 2));
                            winv_aug(icachansind, :) = EEG.icawinv;
                            winv_aug(exclIdx,    :) = W_excl;
                            EEG.icawinv = winv_aug;
                            % Fold sphere into weights and extend to all channels so that
                            % eeg_checkset is satisfied: size(icaweights,2) == size(icasphere,2)
                            %                            == numel(icachansind).
                            % icaweights_aug * eye * data(1:nbchan,:) reproduces icaact
                            % because the excluded-channel columns are zero.
                            weights_aug = zeros(size(EEG.icaweights, 1), EEG.nbchan);
                            weights_aug(:, icachansind) = EEG.icaweights * EEG.icasphere;
                            EEG.icaweights  = weights_aug;
                            EEG.icasphere   = eye(EEG.nbchan);
                            EEG.icachansind = 1:EEG.nbchan;
                        end
                    else
                        % Per-experiment ICA: indices stored directly
                        EEG.icachansind = w.channels;
                        EEG.icasphere   = w.sphere;
                        EEG.icaweights  = w.weights;
                        EEG.icawinv     = w.winverse;
                        EEG.icaact = icaact(EEG.data, EEG.icaweights * EEG.icasphere, mean(EEG.data, 2));
                    end                
            end
            EEG= eeg_checkset(EEG);
        end        
 end


%% Process MFF/Neurostim metadata shared by continuous and epoch datasets.
EEG = addNeurostimMetadata(EEG,key,mffFile,urSrate,false);
[EEG.filepath, name, ext] = fileparts(char(mffFile));
EEG.filename = [name ext];
EEG.etc.neurostim.expt = key;

if EEG.srate ~=urSrate
    for iEvent=1:length(EEG.event)
        EEG.event(iEvent).latency = round(EEG.event(iEvent).latency*(EEG.srate/urSrate));
    end
end
%%
% Check consistency
EEG = eeg_checkset(EEG);

if ~isempty(pv.plg)
    EEG= ephys.egi.eeglabAddEvents(EEG,pv.plg,pv.prm);
end
end
function EEG = epochDataset(key,pv,mffFile)
% Construct an EEGLAB dataset directly from ns.EpochChannel.
epochKey = key;
epochKey.ctag = char(pv.data);
epochKey.etag = char(pv.etag);
if pv.dimension ~= ""
    epochKey.dimension = char(pv.dimension);
end
epochRel = ns.Epoch & epochKey;
nrEpochRows = count(epochRel);
assert(nrEpochRows==1,"Expected exactly one ns.Epoch row for etag=%s; found %d.",pv.etag,nrEpochRows);
epochTime = fetch1(epochRel,'time');
assert(numel(epochTime)==3 && epochTime(3)>=1 && epochTime(3)==round(epochTime(3)),'ns.Epoch.time must be [start stop nrSamples].');
t = linspace(epochTime(1),epochTime(2),epochTime(3));
nrSamples = numel(t);
[channel,trial,signal,onset] = fetchn(ns.EpochChannel & epochRel,'channel','trial','signal','onset');
assert(~isempty(signal),'The selected ns.Epoch contains no EpochChannel rows.');
channel = double(channel(:));
trial = double(trial(:));
if ~iscell(signal), signal = num2cell(signal,2); end
signal = signal(:);
onset = double(onset(:));
[~,order] = sortrows([trial channel],[1 2]);
channel = channel(order); trial = trial(order); signal = signal(order); onset = onset(order);
trialValues = unique(trial,'stable');
channelValues = unique(channel,'stable');
nrTrials = numel(trialValues); nrChannels = numel(channelValues);
data = zeros(nrChannels,nrSamples,nrTrials,'like',signal{1}); seen = false(nrChannels,nrTrials);
for i = 1:numel(signal)
    y = signal{i};
    assert(isvector(y) && numel(y)==nrSamples,'EpochChannel signal for trial %d/channel %d has %d samples; expected %d.',trial(i),channel(i),numel(y),nrSamples);
    if isempty(data), data = zeros(nrChannels,nrSamples,nrTrials,'like',y); end
    ti = find(trialValues==trial(i),1); ci = find(channelValues==channel(i),1);
    assert(~seen(ci,ti),'Duplicate EpochChannel row for trial %d/channel %d.',trial(i),channel(i));
    data(ci,:,ti) = reshape(y,1,[]); seen(ci,ti) = true;
end
assert(all(seen,'all'),'EpochChannel does not contain a complete trial-by-channel grid.');
cRel = ns.C & key & struct('ctag',char(pv.data));
C = fetch(ns.CChannel & cRel,'channel','channelinfo');
cChannels = double([C.channel]);
[isPresent,channelIndex] = ismember(channelValues,cChannels);
assert(all(isPresent),'Epoch channels are missing from the corresponding ns.CChannel relation.');
chanlocs = [C(channelIndex).channelinfo];
EEG = eeg_emptyset();
EEG.setname = sprintf('%s@%sT%s_%s',key.subject,key.session_date,key.starttime,pv.etag);
EEG.srate = round(cRel.samplingRate); EEG.nbchan = nrChannels; EEG.pnts = nrSamples; EEG.trials = nrTrials;
EEG.data = data; EEG.chanlocs = chanlocs; EEG.urchanlocs = chanlocs;
EEG.xmin = t(1); EEG.xmax = t(end); EEG.times = 1000*t;
EEG.event = struct('type',{},'latency',{},'epoch',{},'trial',{});
align = fetch1(ns.EpochParm & epochRel,'align');
assert(isfield(align,'event'),'ns.EpochParm.align must contain an event field.');
alignEvent = char(string(align.event));
alignLatency = -1000*t(1); % Alignment event position within the epoch, in ms.
EEG.epoch = repmat(struct('trial',[],'event',[],'eventtype',alignEvent,'eventlatency',alignLatency,'condition',[],'onset',[]),1,nrTrials);
for i = 1:nrTrials
    EEG.epoch(i).trial = trialValues(i); EEG.epoch(i).onset = onset(find(trial==trialValues(i),1));
end
try
    [dimTrial,condition] = fetchn(ns.DimensionTrial & epochRel,'trial','name');
    for i = 1:nrTrials
        j = find(double(dimTrial)==trialValues(i),1);
        if ~isempty(j), EEG.epoch(i).condition = condition{j}; end
    end
catch
    % Condition metadata is supplementary.
end
EEG.etc.neurostim.expt = key;
EEG.etc.neurostim.epoch = struct('etag',char(pv.etag),'dimension',char(pv.dimension),'time',t,'trials',trialValues,'alignEvent',alignEvent);
EEG = createEpochAlignmentEvents(EEG,alignEvent,trialValues);
[EEG.filepath,name,ext] = fileparts(char(mffFile)); EEG.filename = [name ext];
EEG = eeg_checkset(EEG,'makeur');
end

function EEG = addNeurostimMetadata(EEG,key,mffFile,urSrate,isEpoch)
% Add MFF events, Neurostim trial metadata, clock mapping, and plugin data.
if isEpoch
    [~,begTime] = mff_importinfo(mffFile);
    mffHeader = ft_read_header(mffFile,'headerformat','egi_mff_v1');
    urSrate = mffHeader.Fs;
    EEG.event = mff_importevents(mffFile,begTime,urSrate,0);
end
nrEvts = numel(EEG.event);
eventCode = string({EEG.event.code});
brec = EEG.event(strcmpi('BREC',eventCode));
assert(~isempty(brec),"No BREC event found in " + mffFile + ". Cannot match this EGI file to Neurostim");
EEG.event(strcmpi('BREC',eventCode)).mffkey_TRIA = '1';
[fldr,nsFile,~] = fileparts(file(ns.Experiment & key));
if ~contains(brec.mffkey_FLNM,nsFile)
    jsonFile = fullfile(fldr,nsFile + ".json");
    if exist(jsonFile,"file")
        json = readJson(jsonFile);
        originalFilename = fliplr(extractBefore(fliplr(brec.mffkey_FLNM),'\'));
        ok = contains(json.provenance,originalFilename);
    else
        ok = false;
    end
    assert(ok,sprintf('The MFF file (%s) was created by a different Neurostim file (%s)',brec.mffkey_FLNM,nsFile));
end
isBeginTrial = strcmpi(eventCode,'BTRL');
trial = nan(nrEvts,1);
trial(isBeginTrial) = cellfun(@str2num,{EEG.event(isBeginTrial).mffkey_TRIA});
trial(1) = 1;
trial = fillmissing(trial,"previous");
trial = num2cell(trial);
[EEG.event.trial] = deal(trial{:});
eventEgiTime = ([EEG.event.latency]-1)/urSrate;
prms = get(ns.Experiment & key,{'cic','egi'});
trialStartTimeNeurostim = prms.cic.trial.clocktime(2:end);
trialStartTimeNeurostim = trialStartTimeNeurostim(:)';
trialStartTimeEgi = eventEgiTime(isBeginTrial);
assert(numel(trialStartTimeEgi)==numel(trialStartTimeNeurostim),'Number of trials mismatched in EGI and NS');
if ~isEpoch
    slack = 1;
    startEEGTime = trialStartTimeEgi(1)-slack;
    if numel(trialStartTimeNeurostim)>1
        stopEEGTime = trialStartTimeEgi(end)+median(diff(trialStartTimeEgi))+slack;
        EEG = pop_select(EEG,'time',[startEEGTime stopEEGTime]);
        trialStartTimeEgi = trialStartTimeEgi-startEEGTime;
        EEG.etc.neurostim.clockParms = polyfit(trialStartTimeEgi,trialStartTimeNeurostim,1);
    else
        EEG = pop_select(EEG,'time',[startEEGTime EEG.xmax]);
        trialStartTimeEgi = trialStartTimeEgi-startEEGTime;
        EEG.etc.neurostim.clockParms = [1000 trialStartTimeNeurostim-trialStartTimeEgi*1000];
    end
    % Hack; pop_select can add a boundary event which messes up the prep
    % pipeline later. Delete it.
    if strcmpi(EEG.event(1).type,'boundary'); EEG.event(1)= [];end
else
    EEG.etc.neurostim.clockParms = polyfit(trialStartTimeEgi,trialStartTimeNeurostim,1);
end
% pop_select may have removed some events; reconstruct.
eventEgiTime = ([EEG.event.latency]-1)/urSrate;
eventNsTime = polyval(EEG.etc.neurostim.clockParms,eventEgiTime);
eventTrial = [EEG.event.trial];
eventTrialTime = eventNsTime - trialStartTimeNeurostim(eventTrial);
eventCode = string({EEG.event.code});

if isEpoch
    retainedTrials = [EEG.epoch.trial];
    epochStart = EEG.xmin*1000;
    epochStop = EEG.xmax*1000;
    epochOnset = [EEG.epoch.onset]*1000;
    keep = false(1,numel(EEG.event));
    for i = 1:numel(EEG.event)
        j = find(retainedTrials==eventTrial(i),1);
        if ~isempty(j)
            relativeTime = eventTrialTime(i)-epochOnset(j);
            keep(i) = relativeTime>=epochStart && relativeTime<=epochStop;
            if keep(i)
                EEG.event(i).latency = 1+(relativeTime-epochStart)/1000*EEG.srate;
                EEG.event(i).epoch = j;
            end
        end
    end
    EEG.event = EEG.event(keep);
    eventCode = eventCode(keep);
    eventTrial = eventTrial(keep);
    eventNsTime = eventNsTime(keep);
    eventTrialTime = eventTrialTime(keep);
    for iEpoch = 1:EEG.trials
        EEG.epoch(iEpoch).event = find([EEG.event.epoch] == iEpoch);
    end
end
uNames = unique(eventCode);
prmTpl = struct('property_name','','property_time',[],'property_nstime',[],'property_trial',[],'property_value',[],'property_type','Event');
prmTpl = repmat(prmTpl,[numel(uNames) 1]);
for iName = 1:numel(uNames)
    prmTpl(iName).property_name = char(uNames(iName));
    stay = eventCode == uNames(iName);
    prmTpl(iName).property_time = eventTrialTime(stay);
    prmTpl(iName).property_nstime = eventNsTime(stay);
    prmTpl(iName).property_trial = eventTrial(stay);
    prmTpl(iName).property_value = EEG.event(stay);
end
EEG.etc.neurostim.pluginparameter = prmTpl;
if isEpoch
    % Guarantee one explicit alignment event for every retained epoch.
    alignEventName = EEG.etc.neurostim.epoch.alignEvent;
    alignEvents = repmat(struct('type',alignEventName,'code',alignEventName, ...
        'latency',1+(-EEG.xmin)*EEG.srate,'epoch',0,'trial',0),1,EEG.trials);
    for iEpoch = 1:EEG.trials
        alignEvents(iEpoch).epoch = iEpoch;
        alignEvents(iEpoch).trial = EEG.epoch(iEpoch).trial;
    end
    if isempty(EEG.event)
        EEG.event = alignEvents;
    else
        allFields = union(fieldnames(EEG.event),fieldnames(alignEvents));
        for iField = 1:numel(allFields)
            fieldName = allFields{iField};
            if ~isfield(EEG.event,fieldName)
                [EEG.event.(fieldName)] = deal([]);
            end
            if ~isfield(alignEvents,fieldName)
                [alignEvents.(fieldName)] = deal([]);
            end
        end
        EEG.event = [EEG.event alignEvents];
    end
    for iEpoch = 1:EEG.trials
        EEG.epoch(iEpoch).event = find([EEG.event.epoch] == iEpoch);
    end
end
end
function EEG = ensureEpochAlignmentEvents(EEG)
% Ensure every retained epoch has one EEGLAB alignment event.
if isempty(EEG.event)
    alignEventName = EEG.etc.neurostim.epoch.alignEvent;
    EEG.event = repmat(struct('type',alignEventName,'code',alignEventName, ...
        'latency',1+(-EEG.xmin)*EEG.srate,'epoch',0,'trial',0),1,EEG.trials);
    for iEpoch = 1:EEG.trials
        EEG.event(iEpoch).epoch = iEpoch;
        EEG.event(iEpoch).trial = EEG.epoch(iEpoch).trial;
    end
end
for iEpoch = 1:EEG.trials
    EEG.epoch(iEpoch).event = find([EEG.event.epoch] == iEpoch);
end
end
function EEG = createEpochAlignmentEvents(EEG,alignEvent,trialValues)
% Create one alignment event per epoch from ns.Epoch metadata.
latency = 1+(-EEG.xmin)*EEG.srate;
EEG.event = repmat(struct('type',alignEvent,'code',alignEvent, ...
    'latency',latency,'epoch',0,'trial',0),1,EEG.trials);
for iEpoch = 1:EEG.trials
    EEG.event(iEpoch).epoch = iEpoch;
    EEG.event(iEpoch).trial = trialValues(iEpoch);
    EEG.epoch(iEpoch).event = iEpoch;
end
end