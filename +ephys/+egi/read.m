function  [signal,time,channelInfo,recordingInfo] = read(key,parms)
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

%% Read events
BEGINTIME = 0;
EGIEVENT_SAMPLINGRATE =1000;
evts = mff_importevents(mffFilename, BEGINTIME, EGIEVENT_SAMPLINGRATE); % from time =0 with 1Khz sampling rate
code = {evts.code};
brec = evts(strcmpi('BREC',code)); % neurostim sends this BREC event
assert(~isempty(brec),"No BREC event found in " +  mffFilename + ". Cannot match this EGI file to Neurostim");
evts(strcmpi('BREC',code)).mffkey_TRIA= '1'; % Force it to be in TRIAL 1 (not defined)
[~,nsFile,~] = fileparts(file(ns.Experiment &key));
 % Check that this MFF file was created by the current neurostim file.
assert(contains(brec.mffkey_FLNM,nsFile),sprintf('The MFF file (%s) was created by a different Neurostim file (%s)',brec.mffkey_FLNM,nsFile));

%% Preprocess the events to get trial and time.
nrEvts = numel(evts);
isBeginTrial =strcmpi({evts.code},'BTRL');
trial = nan(nrEvts,1);
trial(isBeginTrial) = cellfun(@str2num,{evts(isBeginTrial).mffkey_TRIA});
trial = fillmissing(trial,"previous");
trial =num2cell(trial);
evts(1).trial=NaN;
[evts.trial]=deal(trial{:});

firstTrialStartTimeEgi = evts(strcmpi(code,'BTRL') & [evts.trial]==1).begintime;
firstTrialStartTimeEgi  = datetime(firstTrialStartTimeEgi(:,1:26),'InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSSSSS');

% trialTime = str2num(fillmissing(char(tcpEvents.mffkey_TTIM),'constant',['-Inf' fillInSpaces])); %#ok<ST2NM> % We place all events without neurostim trial time at -inf trial time.
% data    = egiTime(isTcp); % Seconds since start of file.
% Determine the time of all events (on the EGI clock)
egiTime = char(evts.begintime);
egiTime = datetime(egiTime(:,1:26),'InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
egiTime =num2cell(egiTime);
evts(1).time = NaT;
[evts.time] = deal(egiTime{:});

source= {evts.sourcedevice};
% Split DIN (diode) and TCP (Neurostim) events
isDin = strncmpi(code,'DIN',3);
isTcp = contains(source,'Multi-Port ECI');

%% Read the properties of cic and the egi plugin for this experiment
prms  = get(ns.Experiment & key,{'cic','egi'});
trialStartTimeNeurostim  = prms.cic.trialNsTime(2:end);% 
trialStartTimeEgi = [evts(strcmpi(code,'BTRL')).time];
% These should match
assert(numel(trialStartTimeEgi)==numel(trialStartTimeNeurostim),'Number of trials mismatched in EGI and NS');

if numel(trialStartTimeNeurostim)>10
    clockParms = polyfit((trialStartTimeNeurostim-trialStartTimeNeurostim(1)),1000*seconds(trialStartTimeEgi-trialStartTimeEgi(1)),1);
    fprintf(['Average Clock drift is ' num2str((clockParms(1)-1)) ' ms/ms and the offset is ' num2str(clockParms(2)) ' ms \n' ]);
else
    fprintf('Not enough trials to check clock drift');    
end

%% Read the raw signals from the mff.
fprintf("Start reading from " +  mffFilename + "...")
EEG = mff_import(char(mffFilename)); 
fprintf('Done.\n');
signal =EEG.data';
time = trialStartTimeNeurostim(1)/1000+(0:size(signal,1)-1)/EEG.srate;
%% Preprocess if requested
if isfield(parms,'downsample') || isfield(parms,'designfilt')    
    [signal,time] = ns.CFilter(signal,time,parms);
end
% Regular sampling
time = [ 1000*time(1) 1000*time(end) numel(time)];
signal = single(signal);
channelInfo  = EEG.chanlocs';
channelInfo(1).nr =NaN;
nr = cellfun(@(x) str2double(extractAfter(x,'E')),{channelInfo.labels});
nr = num2cell(nr);
[channelInfo.nr] =deal(nr{:}); 
recordingInfo = mergestruct(EEG.chaninfo,EEG.etc);
recordingInfo.ref = EEG.ref;

%% TODO: Add evts to egi plugin?

end

