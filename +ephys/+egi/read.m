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
% ToDo:
% Add TCP Events (BREC,EREC, BTRL, ETRL, as logged by EGI) and Diode DIN events 
% (As logged by EGI) to the egi plugin.
%  BK - Jan 2025

%% Fetch the file to read (ns.C has already checked that it exists)
mffFilename = strrep(fullfile(folder(ns.Experiment &key),key.filename),'\','/'); % Avoid fprintf errors

%% Read the raw signals from the mff.
fprintf("Start reading from " +  mffFilename + "...")
MFF = mff_import(char(mffFilename)); 
[nrChannels,nrSamples] = size(MFF.data);
fprintf('Done.\n');

%% Read events
code = {MFF.event.code};
brec = MFF.event(strcmpi('BREC',code)); % neurostim sends this BREC event
assert(~isempty(brec),"No BREC event found in " +  mffFilename + ". Cannot match this EGI file to Neurostim");
MFF.event(strcmpi('BREC',code)).mffkey_TRIA= '1'; % Force it to be in TRIAL 1 (not defined)
[~,nsFile,~] = fileparts(file(ns.Experiment &key));
 % Check that this MFF file was created by the current neurostim file.
assert(contains(brec.mffkey_FLNM,nsFile),sprintf('The MFF file (%s) was created by a different Neurostim file (%s)',brec.mffkey_FLNM,nsFile));

%% Preprocess the events to get trial and time.
nrEvts = numel(MFF.event);
isBeginTrial =strcmpi({MFF.event.code},'BTRL');
trial = nan(nrEvts,1);
trial(isBeginTrial) = cellfun(@str2num,{MFF.event(isBeginTrial).mffkey_TRIA});
trial = fillmissing(trial,"previous");
trial =num2cell(trial);
MFF.event(1).trial=NaN;
[MFF.event.trial]=deal(trial{:});

% Determine the time of all events (on the EGI clock)
egiTime = char(MFF.event.begintime);
egiTime = datetime(egiTime(:,1:26),'InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
egiTime =num2cell(egiTime);
MFF.event(1).time = NaT;
[MFF.event.time] = deal(egiTime{:});
source= {MFF.event.sourcedevice};
% Split DIN (diode) and TCP (Neurostim) events
isDin = strncmpi(code,'DIN',3);
isTcp = contains(source,'Multi-Port ECI');

%% Read the properties of cic and the egi plugin for this experiment
% Synchronize clocks. 
% Using NTPSync results in pretty much perfectly aligned clocks (no
% drift,little offset). But we check anyway using the Begin Trial (bTRL)
% events.

prms  = get(ns.Experiment & key,{'cic','egi'});
trialStartTimeNeurostim  = prms.cic.trialNsTime(2:end);% 
trialStartTimeEgi = [MFF.event(strcmpi(code,'BTRL')).time];
% The number of trials should match
assert(numel(trialStartTimeEgi)==numel(trialStartTimeNeurostim),'Number of trials mismatched in EGI and NS');
% Determine clock drift, and the offset between the first trial start event
% in neurostim and in EGI.
brecTimeEgi  = datetime(brec.begintime(1:26),'InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
if numel(trialStartTimeNeurostim)>1
    clockParms = polyfit(seconds(trialStartTimeEgi-brecTimeEgi),trialStartTimeNeurostim,1);
    fprintf(['Average Clock drift is ' num2str((clockParms(1)-1000)) ' ms/s and the offset is ' num2str(clockParms(2)) ' ms \n' ]);
else
    fprintf('Single trial; assuming zero clock drift.\n');     
    % Assume that the trialStartTimeEgi and trialStartTimeNeurostim refer
    % to the same time and that there is no clock drift.
    clockParms = [1000 trialStartTimeNeurostim]; % Assming zero drift, zero offset
end

% EGI samples events regularly, starting from the BREC event 
egiTime = (0:nrSamples-1)/MFF.srate;
neurostimTime = polyval(clockParms,egiTime);
% keep only from the start of the first until the end of the last trial.
stay = neurostimTime >= trialStartTimeNeurostim(1) & neurostimTime <= prms.cic.trialStopTimeNsTime(end);
signal =double(MFF.data(:,stay)');
neurostimTime = neurostimTime(stay);

%% Preprocess if requested
% Apply filtering - time in seconds to allow Hz units for filters
parms.layout = MFF.etc.layout;
[signal,neurostimTime] = ns.CFilter(signal,neurostimTime/1000,parms);
%% Package output
% Regular sampling - stored in ms
neurostimTime = [1000*neurostimTime(1) 1000*neurostimTime(end) numel(neurostimTime)];
signal = single(signal); % Save space -will be converted back to double on read.
channelInfo  = MFF.chanlocs';
channelInfo(1).nr =NaN;
nr = cellfun(@(x) str2double(extractAfter(x,'E')),{channelInfo.labels});
nr = num2cell(nr);
[channelInfo.nr] =deal(nr{:}); 
recordingInfo = mergestruct(MFF.chaninfo,MFF.etc);
recordingInfo.ref = MFF.ref;

%% TODO: Add evts to egi plugin?

end

