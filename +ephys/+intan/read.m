function [signal,time,channelInfo,recordingInfo] = read(key,parms)
% Read from Intan data files, preprocess them according to the parameters
% passed in the parms struct.
% 
% .intan - parameters that determine which channels to read. The neural
% data are called .amplifier. For instance, parms.intan.amplifier = [1 2] will read
% channels 1 and 2. and parms.intan.digIn =NaN will read all digital
% inputs.
% OPTIONAL parms struct members
% .downsample ; The frequency used for downsampling (uses decimate).
% .designfilt ; A struct array of parameters that are passed to designfilt 
%                   and then used to filter the signal with filtfilt.They
%                   will be applied successively (and after downsample; the 
%                   user should not specify the sampling frequency for
%                   designfilt; that is read from the datafile and adjusted
%                   in case of downsampling).
%
% BK - Oct 2023
arguments
    key % The keysource of the Preprocessed table (Experiment and PrepParm tuple)
    parms (1,1) struct  % The preprocessing parameters
end

%% Fetch the file to read (ns.C has already checked that it exists)
filename = fullfile(folder(ns.Experiment &key),fetch1(ns.File &key,'filename'));
[~,~,ext] =fileparts(filename);
%% Read using a script
fprintf('Reading Intan data file %s. ', filename)
switch upper(ext)
    case '.RHS'
        % Call the modified Intan script
        pv  = namedargs2cell(parms.intan);
        [hdr,data]= ephys.intan.read_Intan_RHS2000_file(filename,pv{:},'digIn',NaN);% Read all digIn for the aligment code below
    case '.RHD'
        error('RHD Not implemented yet.')
    otherwise
        error('Unknown Intan data file extension %s .',ext);
end
fprintf('Done in %d seconds.\n ',round(toc))

groups = fieldnames(parms.intan); % amplifier,digIn,digOut, etc.
signal = [];
channelInfo = [];
for i=1:numel(groups)
    if ~isnan(parms.intan.(groups{i})) % Nan means all - no mismatch checking
        nrChannels = numel(parms.intan.(groups{i}));
        assert(size(data.(groups{i}),2)==nrChannels,'Mismatch in requested channels and read channels for %s.',groups{i});
    end
    signal = [signal data.(groups{i})];%#ok<AGROW>
    channelInfo = [channelInfo hdr.(groups{i})]; %#ok<AGROW> -  Per channel filter info.;
end    
[nrSamples,nrChannels]= size(signal);
ch = num2cell(1:nrChannels);
[channelInfo.nr] = deal(ch{:}); % Number them in the order they were concatenated here (if reading different groups together (unlikely ?) , the groups will be concatenated alphabetically)

%% Determine trial start events to  fill the trial mapper and convert Intan time to NS time.
% The Neurostim ripple plugin sends a TTL to the ''trial' digital input channel of the
% intan box and generates a trialStart event; we use this to synchronize
% the clocks.
prms  = get(ns.Experiment & key,{'cic','intanrhx'});
% Extract the time in seconds, on the neurostim clock, when the trialBit was set high.
trialHigh =find([true;diff(prms.intanrhx.trialStartTrial)>0])+1; 
trialStartTimeNeurostim  = prms.intanrhx.trialStartNsTime(trialHigh)/1000- prms.intanrhx.secondsBeforeRead;
nrNsTrials = numel(trialStartTimeNeurostim);
% The neurostim intanrhx plugin names one digIn channel 'trial'
trialChannel = strcmpi({hdr.digIn.custom_channel_name},'trial');
trialStartIx = find([false; diff(data.digIn(:,trialChannel))>0]);  % Down/Up transition
trialStartTimeIntan = data.time(trialStartIx);
nrIntanTrials = numel(trialStartTimeIntan);
if nrIntanTrials==0
    fprintf('No trials in the intan file.. skipping. You may want to remove this file from ns.File\n')
    return
else
    assert(nrNsTrials==nrIntanTrials,'The number of trials in Neurostim (%d) does not match the number recorded by Intan (%d)',nrNsTrials,nrIntanTrials)
end

clockParms =  polyfit(trialStartTimeIntan,trialStartTimeNeurostim,1); % Fit a line to translate intan time to nsTime
resid = polyval(clockParms,trialStartTimeIntan)-trialStartTimeNeurostim;
slope = clockParms(1);
fprintf('Clock residuals: %3.3f ms +/- %3.3f ms, drift %3.3f ms/ms\n',mean(resid),std(resid),slope-1);
if (abs(slope-1)>0.5)
    error('Neurostim-Intan Clock skew larger than 0.5 ms/ms detected');
end
% Convert intan sample time to neurostim time (in seconds)
time  = polyval(clockParms,data.time);


%% Preprocess
if isfield(parms,'downsample') || isfield(parms,'designfilt')
    [signal,time] = ns.filterC(signal,time,parms);
end

recordingInfo = mergestruct(hdr.frequency,hdr.stim);
% Regular sampling so reduce time representation and change to ms.
time = [1000*time(1) 1000*time(end) nrSamples];
% Reduce storage (ns.C.align converts back to double
signal  = single(signal);

end



