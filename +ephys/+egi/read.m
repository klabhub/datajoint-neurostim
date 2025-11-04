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
% parms.eeglab.filt      - Fitering with pop_eegfiltnew. Must be a struct
%                           with fields that match eegfitnew input
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

assert(exist("eeglab.m","file"),'egi.read needs EEGLab. Install it with the MFFMatlabIO and the PREP Pipeline extensions/plugins from https://github.com/sccn/eeglab ')
%% Fetch the file to read (ns.C has already checked that it exists)
exp_tpl = ns.Experiment & key;
mffFilename = strrep(fullfile(folder(exp_tpl),key.filename),'\','/'); % Avoid fprintf errors

%% Read the raw signals from the mff.
fprintf("Using eeglab plugin to read from " +  mffFilename + "...")
% Use the EEG Lab plugin to read the file, with all events.
eegLabSave = 0 ; % Don't save in eeglab
correctEvents = 0; %  Don't correct events with UTF chars/
EEG = pop_mffimport(char(mffFilename),{},eegLabSave,correctEvents);
[nrChannels,nrSamples] = size(EEG.data);
fprintf('Done.\n');

%% Process BREC event to ensure the neurostim and mff file correspond to the same experiment
nrEvts = numel(EEG.event);
eventCode = {EEG.event.code};
brec = EEG.event(strcmpi('BREC',eventCode)); % neurostim sends this BREC event
assert(~isempty(brec),"No BREC event found in " +  mffFilename + ". Cannot match this EGI file to Neurostim");
EEG.event(strcmpi('BREC',eventCode)).mffkey_TRIA= '1'; % Force it to be in TRIAL 1 (not defined)
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
isBeginTrial =strcmpi(eventCode,'BTRL');
trial = nan(nrEvts,1);
trial(isBeginTrial) = cellfun(@str2num,{EEG.event(isBeginTrial).mffkey_TRIA});
trial(1) =1; % Group events before trial 1 with trial 1
trial = fillmissing(trial,"previous");
trial =num2cell(trial);
[EEG.event.trial]=deal(trial{:});
% Determine the time of all events (on the EGI clock)
% Use the .latency field of the MFF.event, not the begintime (which can be
% offset by a few hundred ms).
eventEgiTime = ([EEG.event.latency]-1)/EEG.srate; % This is seconds since the start of the data in EGI
eventEgiTime  =num2cell(eventEgiTime );
[EEG.event.egitime] = deal(eventEgiTime {:});


%% Read the properties of cic and the egi plugin for this experiment
% Synchronize clocks.
% Using NTPSync results in pretty much perfectly aligned clocks (no
% drift). But we check anyway using the Begin Trial (bTRL) events.
prms  = get(exp_tpl,{'cic','egi'});
trialStartTimeNeurostim  = prms.cic.trialNsTime(2:end);%
trialStartTimeEgi = [EEG.event(strcmpi(eventCode,'BTRL')).egitime];
% The number of trials should match
assert(numel(trialStartTimeEgi)==numel(trialStartTimeNeurostim),'Number of trials mismatched in EGI and NS');
% Determine clock drift, and the offset between the first trial start event
% in neurostim and in EGI.
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
egiSampleTime = (0:nrSamples-1)/EEG.srate;
neurostimTime = polyval(clockParms,egiSampleTime);
% Now we have a time axis for the EGI data in Neurostim time
% Keep only from the start of the first until the end of the last trial.
stay = neurostimTime >= trialStartTimeNeurostim(1) & neurostimTime <= prms.cic.trialStopTimeNsTime(end);
channels = 1:nrChannels;

%% Preprocess with EEGLab 
% Fields of the parms.eeglab struct define the operations; they are
% executed in the sequence that they are defined in the parms.eeglab
% struct.
if isfield(parms,'eeglab')
    fn = fieldnames(parms.eeglab);
    for f= 1:numel(fn)
        switch fn{f}
            case 'zapline'
                if isstruct(parms.eeglab.zapline)
                    zapParms = namedargs2cell(parms.eeglab.zapline);
                elseif islogical(parms.eeglab.zapling) && parms.eeglab.zapling
                    zapParms = {};
                end
                EEG = pop_zapline_plus(EEG, zapParms{:});
            case 'prep'
                % Use the PREP pipelin eegLab plugin for preprocessing
                % All unspecified values will be taken from the PREP defaults, except the
                % file paths (which we set to match the MFF file).
                if isstruct(parms.eeglab.prep)
                    % The user specified parameters that differ from the PREP defaults
                    % using a structure with structure fields.
                elseif islogical(parms.eeglab.prep) && parms.eeglab.prep
                    % User had parms.prep =true
                    % Use all PREP defaults
                    %{
        % Not needed documentation only
        parms.prep.general.errorMsgs  = 'verbose';
        parms.prep.boundary.ignoreBoundaryEvents = false;
        parms.prep.resample.resampleOff = true;
        parms.prep.detrend.detrendChannels = 1:nrChannels;
        parms.prep.detrend.detrendType = 'high pass';
        parms.prep.detrend.detrendCutOff = 1;   % Hz
        parms.prep.detrend.detrendStepSize = 0.02; % Seconds
        %parms.prep.globaltrend =- not used.
        parms.prep.linenoise.lineNoiseMethod = 'clean';
        parms.prep.linenoise.lineNoiseChannels = 1:nrChannels;
        parms.prep.linenoise.Fs = EEG.srate;
        parms.prep.linenoise.lineFrequencies = 60:60:round(EEG.srate/2);
        parms.prep.linenoise.p  = 0.01;    
        parms.prep.linenoise.fScanBandwith = 2;
        parms.prep.linenoise.taperBandwith = 2;
        parms.prep.linenoise.taperWindowSize = 4;
        parms.prep.linenoise.taperWindowStep =1;
        parms.prep.linenoise.tau =100;
        parms.prep.linenoise.pad =0 ;
        parms.prep.linenoise.fPassBand = [0 EEG.srate/2];
        parms.prep.linenoise.maximumIterations = 10;
        parms.prep.reference.srate = EEG.srate;
        parms.prep.reference.samples = nrSamples;
        parms.prep.reference.robustDeviationThreshold = 5;
        parms.prep.reference.highFrequencyNoiseThreshold = 5;
        parms.prep.reference.correlationWindowSeconds =1;
        parms.prep.reference.correlationThreshold = 0.4;
        parms.prep.reference.badTimeThreshold = 0.01;
        parms.prep.reference.ransacOff = false;
        parms.prep.reference.ransacSampleSize = 50;
        parms.prep.reference.ransacChannelFraction = 0.25;
        parms.prep.reference.ransacCorrelationThreshold = 0.75;
        parms.prep.reference.ransacUnbrokenTime = 0.4;
        parms.prep.reference.ransacWindowSeconds = 5;
        parms.prep.reference.referenceType = 'robust';
        parms.prep.reference.interpolationOrder = 'post-reference';
        parms.prep.reference.meanEstimateType = 'median';
        parms.prep.reference.referenceChannels = 1:nrChannels;
        parms.prep.reference.evaluationChannels =1:nrChannels;
        % parms.prep.reference.channelLocations = EEG.chanlocs;    %Dont set these
        % parms.prep.reference.channelInformation= EEG.chaninfo;
        parms.prep.reference.maxReferenceIterations =4;
        parms.prep.reference.reportingLevel ='verbose';
        parms.prep.report.reportMode  = 'normal';
        parms.prep.report.summaryFilePath = './summary.pdf';
        parms.prep.report.sessionFilePath = '.report.html';
        parms.prep.report.consoleFID = 1;
        parms.prep.report.publishOn = true;
        parms.prep.report.errorMsgs = 'verbose';
        parms.prep.postprocess.keepFiltered = false;
        parms.prep.postprocess.removeInterpolatedChannels = false;
        parms.prep.postprocess.cleanupReference = false;
                    %}
                    parms.eeglab.prep =struct('report',struct);
                end
                % Only change the reporting file paths if they have not been
                % specified already in the prep parms
                fldr = folder(ns.Session & key);
                if ~isfield(parms.eeglab.prep,'report')
                    parms.eeglab.prep.report = struct('reportMode','normal');
                end
                % Reporting does not work with an absolute file path due to some
                % weirdness in the prep pipeline
                % Temporarily cd
                here= pwd;
                cd(fldr)
                if ~isfield(parms.eeglab.prep.report,'summaryFilePath')
                    parms.eeglab.prep.report.summaryFilePath  = ['.' filesep 'prep_summary.html'];
                end
                if ~isfield(parms.eeglab.prep.report,'sessionFilePath')
                    parms.eeglab.prep.report.sessionFilePath  =  ['.' filesep char(strrep(key.filename,'.mff','_prep.pdf'))];
                end
                try
                    EEG = pop_prepPipeline(EEG, parms.eeglab.prep);
                catch me
                    cd (here)
                    rethrow(me)
                end
                % Keep only the channels that are not marked as stillNoisy
                channels = setdiff(1:nrChannels,EEG.etc.noiseDetection.stillNoisyChannelNumbers)';
                cd (here)
            case 'filt'
                EEG = pop_eegfiltnew(EEG, 'locutoff',parms.eeglab.filtfilt.locutoff,'hicutoff',parms.eeglab.filtfilt.hicutoff,'plotfreqz',0,'usefftfilt',true);
        end
    end
end

signal =EEG.data(channels,stay)';
neurostimTime = neurostimTime(stay);

%% Package output
% Regular sampling - stored in ms
neurostimTime = [neurostimTime(1) neurostimTime(end) numel(neurostimTime)];
signal = single(signal); % Save space on the DJ Server

% RecordingInfo stores information on the session and the preprocessing
recordingInfo.chaninfo = EEG.chaninfo; % All channels - including removed
recordingInfo.etc      = EEG.etc; % Preprocessing results
recordingInfo = makeMymSafe(recordingInfo);

% Channel info stores information per channel (only those that are kept).
channelInfo  = struct2table(EEG.chanlocs(channels));
channelInfo  = addvars(channelInfo,channels,false(height(channelInfo),1),'NewVariableNames',{'nr','interpolated'});
channelInfo.interpolated(ismember(channelInfo.nr,EEG.etc.noiseDetection.interpolatedChannelNumbers)) =true;
channelInfo = makeMymSafe(table2struct(channelInfo));


%% Add evts to egi plugin
% Define the tpls.
plgTpl = fetch(ns.Plugin & 'plugin_name="egi"' & key);

if exists(ns.PluginParameter & plgTpl & 'property_name="BREC"')
    % The events have already been added to the pluginparameter table.
    % (Presumably by some other ns.CParm that also loads from this file).
else
    eventNsTime = polyval(clockParms,[EEG.event.egitime]);
    eventTrialTime = eventNsTime' - trialStartTimeNeurostim([EEG.event.trial]);
    eventCode = string(eventCode);
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

