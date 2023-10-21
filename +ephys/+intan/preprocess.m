function [signal,time,channelInfo,recordingInfo] = preprocess(key,channels,parms)
% Read from Intan data files, preprocess them according to the parameters
% passed in the parms struct.  
%
% LFP , EEG
% parms.designfilt ; these parameters are passed to designfilt and then used to
% filter the signal with filtfilt.
%
% parms.downsample ; used to decimate
%
% NOTE
% Currently only reading/preprocssing amplifier channels (i.e. neural data)
% and using the digIn to align with neurostim. 
% 
% Selective reading of the data files is not trivial due to the way the
% channels are stored. Although skipping all amplifier data could be added
% (for instance to later load/preprocess digIn, or stim data). 
%
arguments
    key % The keysource of the Preprocessed table (Experiment and PrepParm tuple)
    channels (1,:) double {mustBeNonnegative, mustBeInteger} % The list of amplifier channels to preprocess.
    parms (1,1) struct  % The preprocessing parameters
end
qry = ns.File & key & 'extension IN(''.rhs'',''.rhd'')';  % Check for Intan data files
nrFiles = count(qry);
if nrFiles ~=1
    % Zero or more than 1 file
    error('This experiment has %d files. Cannot proceed.');
else
    % Fetch the file to read
    filename = fullfile(folder(ns.Experiment &key),fetch1(qry,'filename'));
end

%% Read file and data
tic
fprintf('Reading Intan data file %s. ', filename)
[hdr, data] = ephys.intan.read(filename,'amplifier',channels,'digIn',NaN);
fprintf('Done in %d seconds.\n ',round(toc))
nrChannels = numel(channels);
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
switch upper(parms.type)    
    case {'EEG','LFP'}       
        % DownSample 
        if isfield(parms,'downsample')
            tic
            fprintf('Downsampling to %.0f Hz (decimate)...',parms.downsample.frequency);    
            R=  round(hdr.frequency.amplifier_sample_rate/parms.downsample.frequency);
            nrSamples = ceil(nrSamples/R);
            tmp = nan(nrSamples,nrChannels);
            for ch = 1:nrChannels
                tmp(:,ch) =  decimate(signal(:,ch),R);
            end
            signal =tmp;
            time = linspace(time(1),time(end),nrSamples)';
            fprintf('Done in %d seconds.\n.',round(toc))
        end
        %% Filtering
        if isfield(parms,'designfilt') && ~isempty(parms.designfilt)
            tic
            fprintf('Applying filter (designfilt)...')       
            d = designfilt(parms.designfilt{:},'SampleRate',hdr.frequency.amplifier_sample_rate);
            signal = filtfilt(d,signal);
            fprintf('Done in %d seconds.\n.',round(toc))
        end        
        
    otherwise 
        error('Unknown Intan preprocessing type %s.',parms.type)
end


channelInfo = hdr.amplifier; % Per channel filter info.
recordingInfo = mergestruct(hdr.frequency,hdr.stim);

end


