function [signal,time,channelInfo,recordingInfo] = readMdaq(key,parms)
% Read from Neurostim .bin files produced by the mdaq plugin.
% Preprocess them according to the parameters passed in the parms struct.
%
% OPTIONAL parms struct members
% .downsample ; The frequency used for downsampling (uses decimate).
% .designfilt ; A struct array of parameters that are passed to designfilt
%                   and then used to filter the signal with filtfilt.They
%                   will be applied successively (and after downsample; the
%                   user should not specify the sampling frequency for
%                   designfilt; that is read from the datafile and adjusted
%                   in case of downsampling).
% .channel ; A cell array of channels that shoudl be read from the mdaq
%               .bin file and stored as continuous data. Channels that
%               don't exist in a given .bin file are skipped. 
%
% .asEvent ; A logical vector indicating for which of the channels an event should be created in the ns.PluginParmaeter table
%           For digital TTL inputs this will mark transitions from low to
%           high or high to low.  For instance for a mdaq channel called
%           'diode', setting asEvent to true will result in the creation of
%           diodeHigh and diodeLow events in the mdaq plugin. These events
%           can then be used to align data ( see ns.Continuous.get())
%
% BK - Oct 2023
arguments
    key % The key (ns.File and ns.ContinuousParm tuple)
    parms (1,1) struct  % The preprocessing parameters
end
if ~iscell(parms.channel)
        parms.channel ={parms.channel};
end

%% Fetch the file to read (ns.Cont has already checked that it exists)
filename = fullfile(folder(ns.Experiment &key),fetch1(ns.File &key,'filename'));


%% Read using a script
tic
load(strrep(filename,'.bin','.mat'),'c');
T=  readBin(c.mdaq,file=filename);
fprintf('Done in %d seconds.\n ',round(toc))

[availableChannels,ix] = intersect(parms.channel,T.Properties.VariableNames);
if numel(availableChannels)~=numel(parms.channel)
        fprintf(2,'Channels %s not recorded in %s; skipped.\n',strjoin(setdiff(parms.channel,availableChannels),'/'),filename); 
        parms.channel= parms.channel(ix);
        parms.asEvent = parms.asEvent(ix);
end

time = seconds(T.nsTime);
signal = T{:,parms.channel};
[nrSamples,nrChannels] = size(signal);
sampleRate = get(c.mdaq.prms.samplerate);
%% Downsample first
if isfield(parms,'downsample')
    tic
    fprintf('Downsampling to %.0f Hz (decimate)...',parms.downsample);
    R=  round(sampleRate/parms.downsample);
    nrSamples = ceil(nrSamples/R);
    tmp = nan(nrSamples,nrChannels);
    for ch = 1:nrChannels
        tmp(:,ch) =  decimate(signal(:,ch),R);
    end
    signal =tmp;
    time = linspace(time(1),time(end),nrSamples)';
    fprintf('Done in %d seconds.\n.',round(toc));
    sampleRate = parms.downsample;
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

%% Create an event for digital ins
if isfield(parms,'asEvent')
    % parms.asEvent is a logical vector that specifies which of the
    % channels is a digital input  and should be represented (also) as an 
    % event in the ns.PluginParameter table.
    digNames = parms.channel(parms.asEvent);
    plgKey  = fetch(ns.Plugin & key & 'plugin_name=''mdaq''');
    digital = signal(:,parms.asEvent)>2.5;
    for ch =1:numel(digNames)        
        % Create an event representing the low to high transition of
        % the digital input
        name =[char(digNames{ch}) 'High'];
        goesHighNsTime = time([false; diff(digital(:,ch))>0])*1000;
        fprintf('Creating %d new %sHigh events in mdaq plugin.\n',numel(goesHighNsTime), digNames{ch})
        [trialTime,trial] = eTime2TrialTime(c.mdaq.prms.vendor,goesHighNsTime);
        nrEvents= numel(goesHighNsTime);    
        if nrEvents>0
            addNew(ns.PluginParameter,plgKey,name,ones(nrEvents,1),'Event',trialTime,trial,goesHighNsTime,true); % Replace existing
        end
        % Create an event representing the high to low transition of
        % the digital input        
        name =[char(digNames{ch}) 'Low'];      
        goesLowNsTime= time([false; diff(digital(:,ch))<0])*1000;
        fprintf('Creating %d new %sLow events in mdaq plugin. \n',numel(goesLowNsTime), digNames{ch})
        [trialTime,trial] = eTime2TrialTime(c.mdaq.prms.vendor,goesLowNsTime);
        nrEvents= numel(goesLowNsTime);
        if nrEvents>0
            addNew(ns.PluginParameter,plgKey,name,ones(nrEvents,1),'Event',trialTime,trial,goesLowNsTime,true);% Replace existing
        end
    end


end

% Information per channel
channelInfo =struct('name',parms.channel,'nr',num2cell(1:nrChannels));
recordingInfo = struct('vendor',c.mdaq.vendor);
% MDaq has regular sampling so reduce time representation
time = [1000*time(1) 1000*time(end) nrSamples];
% Reduce storage (ns.Cont.align converts back to double
signal  = single(signal);

end



