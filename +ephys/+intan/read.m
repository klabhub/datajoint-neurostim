function [signal,time,channelInfo,recordingInfo] = read(key,parms)
% Read from Intan data files, preprocess them according to the parameters
% passed in the parms struct.
% 
% .intan - parameters that determine which channels to read. The neural
% data are called .amplifier. hence parms.intan.amplifier = [1 2] will read
% channels 1 and 2.  See read_Intan_RHS2000_file for the options.
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

%% Check files
qry = ns.File & key;  % Check for Intan data files
nrFiles = count(qry);
if nrFiles ~=1
    % Zero or more than 1 file
    error('This experiment has %d files. Cannot proceed.');
else
    % Fetch the file to read
    filename = fullfile(folder(ns.Experiment &key),fetch1(qry,'filename'));
end

if ~exist(filename,"file")
    error('Intan file %s does not exist',filename)
end
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
nrChannels = numel(parms.intan.amplifier);
assert(size(data.amplifier,2)==nrChannels,'Mismatch in requested channels and read channels.');
signal = data.amplifier;
nrSamples= size(signal,1);

%% Determine trial start events to  fill the trial mapper and convert Intan time to NS time.
% The Neurostim ripple plugin sends a TTL to the ''trial' digital input channel of the
% intan box and generates a trialStart event; we use this to synchronize
% the clocks.
prms  = get(ns.Experiment & key,{'cic','intanrhx'});
% Extract the time in seconds, on the neurostim clock, when the trialBit was set high.
trialStartTimeNeurostim  = prms.intanrhx.trialStartNsTime(find([true;diff(prms.intanrhx.trialStartTrial)>0])+1)/1000;
nrNsTrials = numel(trialStartTimeNeurostim);
% The neurostim intanrhx plugin names one digIn channel 'trial'
trialChannel = strcmpi({hdr.digIn.custom_channel_name},'trial');
trialStartIx = find([false; diff(data.digIn(:,trialChannel))>0]);  % Down/Up transition
trialStartTimeIntan = data.time(trialStartIx);
nrIntanTrials = numel(trialStartTimeIntan);
assert(nrNsTrials==nrIntanTrials,'The number of trials in Neurostim (%d) does not match the number recorded by Intan (%d)',nrNsTrials,nrIntanTrials)

clockParms =  polyfit(trialStartTimeIntan,trialStartTimeNeurostim,1); % Fit a line to translate intan time to nsTime
resid = polyval(clockParms,trialStartTimeIntan)-trialStartTimeNeurostim;
slope = clockParms(1);
fprintf('Clock residuals: %3.3f ms +/- %3.3f ms, drift %3.3f ms/ms\n',mean(resid),std(resid),slope-1);
if (abs(slope-1)>0.5)
    error('Neurostim-Intan Clock skew larger than 0.5 ms/ms detected');
end
% Convert intan sample time to neurostim time
time  = polyval(clockParms,data.time);


%% Preprocess

%% Downsample first
if isfield(parms,'downsample')
    tic
    fprintf('Downsampling to %.0f Hz (decimate)...',parms.downsample);
    R=  round(hdr.frequency.amplifier_sample_rate/parms.downsample);
    nrSamples = ceil(nrSamples/R);
    tmp = nan(nrSamples,nrChannels);
    for ch = 1:nrChannels
        tmp(:,ch) =  decimate(signal(:,ch),R);
    end
    signal =tmp;
    time = linspace(time(1),time(end),nrSamples)';
    fprintf('Done in %d seconds.\n.',round(toc));
    sampleRate = parms.downsample;
else
    sampleRate = hdr.frequency.amplifier_sample_rate;
end
%% Notch, Bandpass,etc.
if isfield(parms,'designfilt')
    fn = fieldnames(parms.designfilt);
    for i=1:numel(fn)
        tic;
        fprintf('Applying filter (designfilt.%s)...',fn{i})
        prms= parms.designfilt.(fn{i});
        d = designfilt(prms{:},'SampleRate',sampleRate);
        signal = filtfilt(d,signal);
        fprintf('Done in %d seconds.\n.',round(toc))
    end
end
channelInfo = hdr.amplifier; % Per channel filter info.
ch = num2cell(1:nrChannels);
[channelInfo.nr] = deal(ch{:}); % Number them in the order they were requested (as parms.amplifier)
recordingInfo = mergestruct(hdr.frequency,hdr.stim);

end



