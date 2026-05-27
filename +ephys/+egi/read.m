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

EEG = ephys.egi.eeglabDataset(key,data="RAW");




%% Preprocess with EEGLab
% Fields of the parms.eeglab struct define the operations; they are
% executed in the sequence that they are defined in the parms.eeglab
% struct.

fn = fieldnames(parms.eeglab);
for f= 1:numel(fn)
    switch fn{f}
        case 'resample'
            if iscell(parms.eeglab.resample)
                % Passed verbatim to pop_resample
                % Must be {freq [Hz] ,fc,df}
                %   fc         - anti-aliasing filter cutoff (pi rad / sample)
                %                {default 0.9}
                %   df         - anti-aliasing filter transition band width (pi rad /
                %                sample) {default 0.2}
                % fc and df are optional
                resampleParms = parms.eeglab.resample;
            elseif isnumeric(parms.eeglab.resample) && isscalar(parms.eeglab.resample)
                % Only frequency specified.
                resampleParms = {parms.eeglab.resample};
            else
                error('parms.eeglab.resample must be either a scalar double (frequency) or a cell array with frequency,fc, and df. see pop_resample');
            end
            EEG = pop_resample(EEG, resampleParms{:});
        case 'zapline'
            if isstruct(parms.eeglab.zapline)
                if ~isfield(parms.eeglab.zapline,'noisefreqs')
                    % Oddly zapline does not look at line (50/60) noise
                    % by default. Force that default here
                    parms.eeglab.zapline.noisefreqs = 'line';
                end
                zapParms = namedargs2cell(parms.eeglab.zapline);
            elseif islogical(parms.eeglab.zapline) && parms.eeglab.zapline
                zapParms = {'noisefreqs','line'};
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
            EEG = pop_select(EEG, 'channel', setdiff(1:EEG.nbchan,EEG.etc.noiseDetection.stillNoisyChannelNumbers)');
            cd (here)
        case 'filt'
            EEG = pop_eegfiltnew(EEG, 'locutoff',parms.eeglab.filt.locutoff,'hicutoff',parms.eeglab.filt.hicutoff,'plotfreqz',0,'usefftfilt',true);
        otherwise
            error('Unknown eeglab preprocessing struct %s  \n',fn{f})
    end
end


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

