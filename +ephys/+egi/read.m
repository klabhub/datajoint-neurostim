function  [signal,neurostimTime,channelInfo,recordingInfo] = read(key,parms)
% Function to read MFF files created by EGI. The reading itself uses the
% mffmatlabio plugin by Arno Delorme installed as a plugin in EEGLab.
% EEGLab must be on the search path.
%
% File reading relies on the convention that all EGI mff files are named
% subject.paradigm_day_time.mff and stored in the same directory as the
% Neurostim output file.
%
% During the experiment, the experimenter types subject.paradigm in the Netstation
% box asking for the subject ID, and Netstation adds day_time. As a
% consequence, the time will not necessarily match the time of the
% Neurostim file. Neurostim, however, sends a _FLNM event with the neurostim
% file name. We use this to match NS and MFF files.
%
% The parms.eeglab struct defines preprocessing operations done with eeglab
% and its plugins. Currently we use
% parms.eeglab.zapline   - Zapline-plus line noise removal
%           Set this to true to use all default parameters, or set it to a
%           struct with fieldnames that match the zapline parameters to use
%           different settings.
% parms.eeglab.prep      - the PREP pipeline
%           Set this to true to use all default parameters, or set it to a
%           struct with fieldnames that match the PREP parameters to use
%           different settings.
% parms.eeglab.filtfilt  - Filtering with pop_eegfiltnew. Must be a struct
%                           with fields that match eegfiltnew input
%                           parameters.
%
% The TCP events in the MFF file  contain timing events sent by Neurostim (see
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
% Nov 25 -  Moved to using EEGLab and its plugins.
arguments
    key (1,1) struct
    parms (1,1) struct
end
assert(exist("eeglab.m","file"),'egi.read needs EEGLab. Install it with the MFFMatlabIO and the PREP Pipeline extensions/plugins from https://github.com/sccn/eeglab ')
assert(isfield(parms,'eeglab'),'Preprocessing of EGI MFF files must use EEGLAB. Define the parms.eeglab struct.');

% Read from the MFF file
EEG = ephys.eeglab.dataset(key,data="RAW");
% Run preprocessing
EEG = ephys.eeglab.preprocess(EEG,parms); 


%% Package output and save full results to file.
% The output in EEG.etc can be huge (the prepline for instance stores
% some of the raw data  in there). Save that to a file and keep only a
% few fields in the database.
mffFilename = fullfile(EEG.filepath,EEG.filename);
etcFile = strrep(mffFilename,'.mff','_etc.mat');
etc = EEG.etc;
w = whos('etc');
fprintf('Saving %.1f MB preprocessing results to %s\n',w.bytes/1e6,etcFile);
save(etcFile,'etc')

% Select a subset of results to store in the database. (the rest can
% in principle be retrievd from the etcFile later).

% RecordingInfo stores information on the session and the preprocessing
recordingInfo.chaninfo = EEG.chaninfo; % All channels - including removed
% Channel info stores information per channel (only those that are kept).
channelInfo  = struct2table(EEG.chanlocs);
channelInfo.nr = cellfun(@str2double,extractAfter(channelInfo.labels,'E')); % Keep original numbering

if isfield(etc,'noiseDetection')
    % Prep pipeline was used
    recordingInfo.etc.noiseDetection.version = etc.noiseDetection.version;
    recordingInfo.etc.noiseDetection.interpolatedChannelNumbers = etc.noiseDetection.interpolatedChannelNumbers;
    recordingInfo.etc.noiseDetection.removedChannelNumbers = etc.noiseDetection.removedChannelNumbers;
    recordingInfo.etc.noiseDetection.stillNoisyChannelNumbers = etc.noiseDetection.stillNoisyChannelNumbers;
    channelInfo  = addvars(channelInfo,false(height(channelInfo),1),'NewVariableNames',{'interpolated'});
    channelInfo.interpolated(ismember(channelInfo.nr,EEG.etc.noiseDetection.interpolatedChannelNumbers)) =true;
end


neurostimTime = polyval(EEG.etc.neurostim.clockParms, EEG.times/1000);
neurostimTime = [neurostimTime(1) neurostimTime(end) numel(neurostimTime)];
signal = single(EEG.data'); % Save space on the DJ Server
recordingInfo = makeMymSafe(recordingInfo);
channelInfo = makeMymSafe(table2struct(channelInfo));


%% Add evts to egi plugin
% Define the tpls.
plgTpl = fetch(ns.Plugin & 'plugin_name="egi"' & key);
if exists(ns.PluginParameter & plgTpl & 'property_name="BREC"')
    % The events have already been added to the pluginparameter table.
    % (Presumably by some other ns.CParm that also loads from this file).
else
    prmTpl= mergestruct(plgTpl,EEG.etc.neurostim.pluginparameter);
    insert(ns.PluginParameter,prmTpl);
end
end

