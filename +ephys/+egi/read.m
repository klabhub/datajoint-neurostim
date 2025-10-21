function  [signal,neurostimTime,channelInfo,recordingInfo] = read(key,parms)
% Function to read MFF files created by EGI. The reading itself uses the
% mffmatlabio toolbox by Arno Delorme  (https://github.com/arnodelorme/mffmatlabio).
% This toolbox must be on the  matlab search path.
%
% This relies on the convention that all EGI mff files are named
% subject.paradigm_day_time.mff and stored in the same directory as the
% Neurostim output file.
%
% During the experiment, the experimenter types subject.paradigm in the Netstation
% box asking for the subject ID, and Netstation adds day_time. As a
% consequence, the time will not necessarily match the time of the
% Neurostim file. Neurostim, however, sends a _FLNM event with the neurostim
% file name. We use this to match NS and MFF files.
%
% The TCP events in the MFF file will contain timing events sent by Neurostim (see
% neurostim.plugins.egi).
%
% Conventions (all files in the same diretory as the neurostim mat file):
%   Main MFF File:  subject.paradigm.day_time.mff
%   GPS Mff File:    subject.gps.day_time.mff
%   Solved coordinates from GPS:   subject.coordinates.xml
%   Impedance checks:   subject.zcheck.day_time.mff
%
%  BK - Jan 2025
%
% Oct 25 - Fixed a bug in time alignment by using MFF.event.latency instead
% of .begintime
%       - Added events to plugin parameter table

%% Fetch the file to read (ns.C has already checked that it exists)
exp_tpl = ns.Experiment & key;
mffFilename = strrep(fullfile(folder(exp_tpl),key.filename),'\','/'); % Avoid fprintf errors

%% Read the raw signals from the mff.
fprintf("Start reading from " +  mffFilename + "...")
MFF = mff_import(char(mffFilename));
[nrChannels,nrSamples] = size(MFF.data); %#ok<ASGLU>
fprintf('Done.\n');

%% Read events
code = {MFF.event.code};
brec = MFF.event(strcmpi('BREC',code)); % neurostim sends this BREC event
assert(~isempty(brec),"No BREC event found in " +  mffFilename + ". Cannot match this EGI file to Neurostim");
MFF.event(strcmpi('BREC',code)).mffkey_TRIA= '1'; % Force it to be in TRIAL 1 (not defined)
[fldr,nsFile,~] = fileparts(file(ns.Experiment &key));
% Check that this MFF file was created by the current neurostim file.
if ~contains(brec.mffkey_FLNM,nsFile)
    % No match. Check if the file was renamed with nsMeta
    jsonFile = fullfile(fldr,nsFile + ".json");
    if exist(jsonFile,"file")
        json = readJson(jsonFile);
        originalFilename = fliplr(extractBefore(fliplr(brec.mffkey_FLNM),'\'));
        ok= contains(json.provenance,originalFilename);  % OK: this was renamed after recording
    else
        ok =false;
    end
    assert(ok ,sprintf('The MFF file (%s) was created by a different Neurostim file (%s)',brec.mffkey_FLNM,nsFile));
end

%% Preprocess the events to get trial and time.
nrEvts = numel(MFF.event);
isBeginTrial =strcmpi({MFF.event.code},'BTRL');
trial = nan(nrEvts,1);
trial(isBeginTrial) = cellfun(@str2num,{MFF.event(isBeginTrial).mffkey_TRIA});
trial(1) =1; % Group events before trial 1 with trial 1
trial = fillmissing(trial,"previous");
trial =num2cell(trial);
[MFF.event.trial]=deal(trial{:});

% Determine the time of all events (on the EGI clock)
% Use the .latency field of the MFF.event, not the begintime (which can be
% offset by a few hundred ms).
eventEgiTime = ([MFF.event.latency]-1)/MFF.srate; % This is seconds since the start of the data in EGI
eventEgiTime  =num2cell(eventEgiTime );
[MFF.event.egitime] = deal(eventEgiTime {:});


%% Read the properties of cic and the egi plugin for this experiment
% Synchronize clocks.
% Using NTPSync results in pretty much perfectly aligned clocks (no
% drift,little offset). But we check anyway using the Begin Trial (bTRL)
% events.
prms  = get(exp_tpl,{'cic','egi'});
trialStartTimeNeurostim  = prms.cic.trialNsTime(2:end);%
trialStartTimeEgi = [MFF.event(strcmpi(code,'BTRL')).egitime];
% The number of trials should match
assert(numel(trialStartTimeEgi)==numel(trialStartTimeNeurostim),'Number of trials mismatched in EGI and NS');
% Determine clock drift, and the offset between the first trial start event
% in neurostim and in EGI.
%brecTimeEgi  = datetime(brec.begintime(1:26),'InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
if numel(trialStartTimeNeurostim)>1
    clockParms = polyfit(trialStartTimeEgi,trialStartTimeNeurostim,1);
    fprintf(['Average Clock drift is ' num2str((clockParms(1)-1000)) ' ms/s and the offset is ' num2str(clockParms(2)) ' ms \n' ]);
else
    fprintf('Single trial; assuming zero clock drift.\n');
    % Assume that the trialStartTimeEgi and trialStartTimeNeurostim refer
    % to the same time and that there is no clock drift. The 1000 slope
    % takes the s from egi to ms used here.
    clockParms = [1000 trialStartTimeNeurostim]; % Assming zero drift, zero offset
end

% EGI samples events regularly:
egiSampleTime = (0:nrSamples-1)/MFF.srate;
neurostimTime = polyval(clockParms,egiSampleTime);
% Now we have a time axis for the EGI data in Neurostim time
% Keep only from the start of the first until the end of the last trial.
stay = neurostimTime >= trialStartTimeNeurostim(1) & neurostimTime <= prms.cic.trialStopTimeNsTime(end);
signal =double(MFF.data(:,stay)');
neurostimTime = neurostimTime(stay);

%% Preprocess if requested

% Layout necessary for certain referencing and interpolation functions
parms.layout = MFF.etc.layout;
parms.layout.ChannelLocations = [[MFF.chanlocs.X]', [MFF.chanlocs.Y]', [MFF.chanlocs.Z]'];

% If assigned badElectrodes file, the electrodes will have been marked
% during preprocessing
badElectrodes = ephys.egi.badElectrodes(exp_tpl, ...
    struct(filename = 'badElectrodes', extension = '.xlsx'));
if isfield(badElectrodes, "channel")
    parms.badElectrodes = badElectrodes.channel;
    signal(parms.badElectrodes, :) = NaN;
else
    parms.badElectrodes = [];
end
% Pass to CFilter (needs time in s for filter definitions), returns time
% in s.
[signal,neurostimTime,prepResult] = ns.CFilter(signal,neurostimTime/1000,parms,exp_tpl);

%% Package output
% Regular sampling - stored in ms
neurostimTime = [1000*neurostimTime(1) 1000*neurostimTime(end) numel(neurostimTime)];
signal = single(signal); % Save space.
channelInfo  = MFF.chanlocs';
channelInfo(1).nr =NaN;
nr = cellfun(@(x) str2double(extractAfter(x,'E')),{channelInfo.labels});
nr = num2cell(nr);
[channelInfo.nr] =deal(nr{:});
recordingInfo = mergestruct(MFF.chaninfo,MFF.etc);
recordingInfo.ref = MFF.ref;
recordingInfo.srate = MFF.srate; % This is the original rate, not the rate of signal now(which could be downsampled)
recordingInfo.layout = parms.layout;
recordingInfo = mergestruct(recordingInfo,prepResult);

%% Add evts to egi plugin
% Define the tpls.
plgTpl = fetch(ns.Plugin & 'plugin_name="egi"' & key);

if exists(ns.PluginParameter & plgTpl & 'property_name="BREC"')
    % The events have already been added to the pluginparameter table.
    % (Presumably by some other ns.CParm that also loads from this file).
else
    eventNsTime = polyval(clockParms,[MFF.event.egitime]);
    eventTrialTime = eventNsTime' - trialStartTimeNeurostim([MFF.event.trial]);
    eventCode = string({MFF.event.code});
    uNames= unique(eventCode);
    nrNames= numel(uNames);

    prmTpl  = mergestruct(plgTpl,struct('property_name','','property_time',[],'property_nstime',[],'property_trial',[],'property_value',[],'property_type','Event'));
    prmTpl = repmat(prmTpl,[nrNames 1]);

    nmCntr =0;
    for nm = uNames
        nmCntr = nmCntr+1;
        prmTpl(nmCntr).property_name = char(nm);

        stay = eventCode == nm;
        prmTpl(nmCntr).property_time = eventTrialTime(stay);
        prmTpl(nmCntr).property_nstime=eventNsTime(stay);
        prmTpl(nmCntr).property_trial=[MFF.event(stay).trial];
        prmTpl(nmCntr).property_value = [MFF.event(stay)];
    end
    insert(ns.PluginParameter,prmTpl);
end
end

