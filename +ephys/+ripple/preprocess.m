function [signal,timeNeurostim,info] = preprocess(key,channels,parms)
% Read from Ripple data files, preprocess them according to the parameters
% passed in the parms struct.  This handles MUAE, LFP, SPIKES, and EEG
% data and relies on the Ripple neuroshare tools to read the nev, nsx etc
% files.
%
% MUAE:
% see multiUnitAvtivityEnvelope for settings, but there is probably little 
%           reason to change those.
% LFP
% parms.designfilt ; these parameters are passed to designfilt and then used to
% filter the signal with filtfilt.
%
arguments
    key % The keysource of the Preprocessed table (Experiment and PrepParm tuple)
    channels (1,:) double % The list of channels to preprocess.
    parms (1,1) struct % The preprocessing parameters
end
import ephys.ripple.*
% Determine  which label to use
switch upper(parms.type)
    case {'MUAE','RAW'}
        label = 'raw';
    case {'LFP'}
        label = 'lfp';
    case {'EEG'}
        label = 'hi-res';
    case {'BREAKOUT'}
end

qry = ns.File & key & 'extension=''.nev''';  % This should always exist for a Ripple recording
nrFiles = count(qry);
if nrFiles ~=1
    % Zero or more than 1 file
    error('This experiment has %d nev-files. Cannot proceed.');
else
    % Fetch the file to read
    filename = fullfile(folder(ns.Experiment &key),fetch1(qry,'filename'));
end

%% Open it with neuroshare
tic
fprintf('Reading header from %s. ', filename)
if ~exist(filename,'file')
    error('ripple:open','File %s does not exist',filename)
else
    [errCode, hFile] =ns_OpenFile(char(filename));
end
if ~strcmpi(errCode,'ns_OK');error('ns_OpenFile failed with %s', errCode);end
fprintf('Done in %d seconds.\n ',round(toc))

entities = [hFile.Entity];

%% Find relevant channels
entityIx = [];
for i=1:numel(entities)
    if strcmpi(entities(i).EntityType,'Analog') 
        thisChannel = extractAfter(entities(i).Label,label);
        if ~isempty(thisChannel) && ismember(str2double(thisChannel),channels)
            entityIx = [entityIx i]; %#ok<AGROW> 
        end
    end
end
nrChannels = numel(entityIx);
fprintf('%d channels from this Array in this file\n',nrChannels)
if nrChannels ==0
    signal = []; timeNeurostim = [];info=[];
    return;
end
assert(nrChannels==numel(channels),'Multiple channel matches?')
nrSamples = unique([entities(entityIx).Count]);
assert(numel(nrSamples)==1,"The code assumes all %s have the same number of samples",label);


%% Preprocess
switch upper(parms.type)
    case 'MUAE'
        for i=1:nrChannels
            [errCode, info(i)] = ns_GetAnalogInfo(hFile, entityIx(i)); %#ok<AGROW>
            if ~strcmpi(errCode,'ns_OK');error('ns_GetAnalogInfo failed with %s', errCode);end
        end
        samplingRate = info(1).SampleRate;% All are  the same for sure.
        rawTime = ((0:nrSamples-1)'/samplingRate);
        nrSamplesInMuae = numel(downsample(rawTime,samplingRate/parms.targetHz));
        signal = nan(nrSamplesInMuae,nrChannels);
        % Read at least one channel at a time, but if maxMemory allows,
        % read more (potentially to speed up processing).
        memoryAvailable = parms.maxMemory; 
        channelChunk = floor(memoryAvailable/(nrSamples*8/1e9))+1;
        for i=1:channelChunk:nrChannels
            tic
            thisChunk = i:min(nrChannels,i+channelChunk-1);
            fprintf('Reading RAW #%d to #%d: channel %d to %d from file ...',min(thisChunk),max(thisChunk),entities(entityIx(min(thisChunk))).ElectrodeID,entities(entityIx(max(thisChunk))).ElectrodeID)
            [errCode, raw] = ns_GetAnalogDataBlock(hFile, entityIx(thisChunk), 1, nrSamples,'scale');
            if ~strcmpi(errCode,'ns_OK');error('ns_GetAnalogDataBlock failed with %s', errCode);end
            fprintf('Converting to MUAE ...')
            [signal(:,thisChunk),timeRipple] = ephys.multiUnitActivityEnvelope(raw,rawTime,"parms",parms);
            fprintf('Done in %d seconds.\n.',round(toc))
        end        
    case {'EEG','LFP'}
        for i=1:nrChannels
            [errCode, info(i)] = ns_GetAnalogInfo(hFile, entityIx(i)); %#ok<AGROW>
            if ~strcmpi(errCode,'ns_OK');error('ns_GetAnalogInfo failed with %s', errCode);end
        end
        samplingRate = info(1).SampleRate;% All are  the same for sure.
        timeRipple = (0:nrSamples-1)/samplingRate;        
        tic
        fprintf('Reading from file...')
        [errCode, signal] = ns_GetAnalogDataBlock(hFile,  entityIx, 1, nrSamples,'scale');
        if ~strcmpi(errCode,'ns_OK');error('ns_GetAnalogDataBlock failed with %s', errCode);end
        fprintf('Done in %d seconds.\n.',round(toc))    
        %% Filtering
        if isfield(parms,'designfilt') && ~isempty(parms.designfilt)
            tic
            fprintf('Applying filter (designfilt)... ')       
            d = designfilt(parms.designfilt{:});
            signal = filtfilt(d,signal);
            fprintf('Done in %d seconds.\n.',round(toc))
        end
        
            

    case 'SPIKES'
        error('Ripple preprocessing for %s has not been implemented yet.',parms.type)
    otherwise 
        error('Unknown Ripple preprocessing type %s.',parms.type)
end


%% Determine trial start events to  fill the trial mapper and convert nip time to NS time.
% The Neurostim ripple plugin sets the digital out of the NIP and generates a trialStart event
% We use this to synchronize the clocks
prms  = get(ns.Experiment & key,{'cic','ripple'});
% The first event is the correct one.. (sic)
% This is the time, on the neurostim clock, when the trialBit was set
% high.
trialStartTimeNeurostim  = prms.ripple.trialStartNsTime(find([true;diff(prms.ripple.trialStartTrial)>0])+1)/1000;
% Find wich bit stored the trialStartEvent and get the time on the NIP
bit = get(ns.Experiment & key,'ripple','prm','trialBit');
eventIx  = find(ismember({entities.EntityType},'Event'));
expression = ['\<SMA\s*' num2str(bit)];
trialBitEntityIx  = find(~cellfun(@isempty,regexp({entities(eventIx).Reason},expression,'match')));
[errCode, trialBitTime,trialBitValue] = ns_GetEventData(hFile, eventIx(trialBitEntityIx), 1:entities(eventIx(trialBitEntityIx)).Count);
if ~strcmpi(errCode,'ns_OK');error('ns_GetEventData failed with %s', errCode);end
% There are various issues with the timing on the NIP (some events can be
% msising, others happen more than once...). This function deals with those
% and then returns the linear parameters that link the two clocks,
clockParms = matchRiplleNeurostim(trialBitTime(:),trialBitValue(:),trialStartTimeNeurostim(:));
% Use this to translate the time of the preprocessed samples (timeRipple)
% to time on the neurostim clock.
timeNeurostim =  polyval(clockParms,timeRipple); % Conver ripple time to nsTime

ns_CloseFile(hFile);

end



function [clockParms] = matchRiplleNeurostim(time,value,neurostimTime,maxSlack)
% Check whether the time between events is the same in ripple and neurostim (and
% correct the start/stop times if Ripple has some extraneous
% bit high events.). The maxSlack (ms) input states which differences in time
% between events is considered too much (and generate an error
% or an attempt to fix the mismatch).
% start = times when the bit went high
% stop  = times when the bit went low
% missing  = trial numbers in which nev trialbits were added
% (i.e. missing from .nev file)
% clockParms = linear fit of the start and startInNs clocks

if nargin <4
    maxSlack = 5; % ms should be enough for events that are not related to the monitor refresh (i.e. trialbit)
end

% Remove zeros at the leading edge.
ix =find(value>0,1);
if isempty(ix)% No nonzero values; happens with UDP loopback.
    return;
end
value(1:(ix-1)) = [];
time(1:(ix-1))=[];
% With UDP loopback enabled, digital output values are stored multiple times.
% Here we detect that and just store the first.
flip = [true; diff(value)~=0];
start = time(value==32767);
stop = time(value==0);
if sum(flip) ~=numel(flip)
    %Seems robust. No need to warn. fprintf('Fixing Trellis UDP Lookback duplicate Up/Down events.\n')
    start = time(flip & value==32767);
    stop  = time(flip & value==0);
end
if length(start) ~= length(stop) || isempty(start)
    if length(stop) == length(start)-1 && stop(end)<start(end)
        fprintf('Last trial stop missing from Ripple file. Estimating.\n');
        stop = [stop;stop(end) + prctile(diff(start),95)];
    else
        fprintf('Start and stop trial information is not complete.\n');
    end
end

[~,~,clockParms] = matchTwoClocks(start,neurostimTime,maxSlack,stop);

end


function [newRippleStart,newRippleStop,parms] = matchTwoClocks(rippleStart,nsStart,slack,rippleStop)
% Match two clocks - returns the ripple times that match the ns
% times.
% rippleStart = start time on clock 1
% nsStart  - start time on clock 2. Clock 2 is assume to be the one to match.
% slack = Slack to allow for matching time between clock ticks on each of
% the clocks.
ADDSLACK =0.5;
nsStart=nsStart(:);
nrNs = numel(nsStart);
rippleStart=rippleStart(:);
rippleStop=rippleStop(:);
nsDuration = diff(nsStart);
rippleDuration = diff(rippleStart);

minNrMatches = 10; % Once you have found this number, do a linear fit.

%% First look for an initial period of matched durations.
% Find the first that matches within slack.
nrEntries = min(numel(rippleDuration),numel(nsDuration)); %
firstTr = find(abs(nsDuration(1:nrEntries)-rippleDuration(1:nrEntries))<slack,1);
firstNonMatch = find(abs(nsDuration(firstTr+1:nrEntries)-rippleDuration(firstTr+1:nrEntries))>slack,1);
if isempty(firstNonMatch)
    % All match up to nrEntries.
    tr = nrEntries;
else
    tr = firstTr+ firstNonMatch; % start looping from here
    while(tr<=nrEntries)
        if abs(nsDuration(tr)-rippleDuration(tr))<slack
            % Good enough match, move on to the next
            tr=tr+1;
        elseif tr<nrEntries && abs(nsDuration(tr)-sum(rippleDuration(tr:tr+1)))<slack
            % Two starst in ripple match one in ns; this
            % probably means the start at tr was an error.
            rippleStop(tr+1)=[];
            rippleStart(tr+1)=[];
            rippleDuration = diff(rippleStart);
            nrEntries = min(numel(rippleDuration),numel(nsDuration)); %
            tr = tr+1;
        else
            nrMatches = tr-firstTr-1; % Current one is not a match so -1
            if nrMatches > minNrMatches
                tr=tr-1; % Move back one
                break; % found enough - use a linear fit
            else
                % increase slack
                slack = slack + ADDSLACK ;
                fprintf('Clock match slack was increased to %.2f\n ', slack);
            end
        end
    end
end
tr = tr +1; % corrct from duration (diff) to trial.

%% Fit a line to sync clocks
parms = polyfit(rippleStart(firstTr:tr),nsStart(firstTr:tr),1);
resid = polyval(parms,rippleStart(firstTr:tr))-nsStart(firstTr:tr);
offset = parms(2);
slope = parms(1);
fprintf('Clock residuals: %3.3f ms +/- %3.3f ms, drift %3.3f ms/ms\n',mean(resid),std(resid),slope-1);
if (abs(slope-1)>0.01)
    error('Ns-Ripple Clock skew larger than 0.01 detected');
end

if tr==nrNs % Found matches all the way to the last trial
    if firstTr==1
        % All matched within slack (extra ripple bits were removed)
    elseif (tr-firstTr+1) == nrNs
        % First few in ripple skipped, but the rest matched
    elseif numel(rippleStart)==nrNs
        % The firstTr did not match within slack, but we
        % probably don't need them anyway
        fprintf('First %d trials have %3.3f ms duration mismatch (included anyway).\n',firstTr-1,mean(rippleDuration(1:firstTr-1)-nsDuration(1:firstTr-1)));
        firstTr = 1;
    end
    % clock skew fine
    newRippleStart = rippleStart(firstTr:tr);
    newRippleStop = rippleStop(firstTr:tr);
else
    % use the fit to match all
    invParms = [1./slope -offset];
    newRippleStart = polyval(invParms,nsStart);
    newRippleStop = newRippleStart + [nsDuration;prctile(nsDuration,99)];
end

end