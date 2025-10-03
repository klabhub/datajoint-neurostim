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
% .analogChannel ; A cell array of channels to read from the mdaq
%               .bin file and stored as continuous data. Channels that
%               don't exist in a given .bin file are skipped.
% .digitalChannel ; A cell array of channels to read with 0/1 entries.
%                   Transitions from low to high or high to low will be added as events.
%                For instance for a mdaq channel called
%            'diode', will result in the creation of
%           diodeHigh and diodeLow events in the mdaq plugin. These events
%           can then be used to align data ( see ns.Experiment.get() and ns.C.align())
% .asEvent Fields in this struct identify analog channels that should be
% thresholded and converted to events. The value of the field is used as
% the threshold, if it is nan, the midpoint between the 5th and 95th
% percentil is used. E.g. .asEvent.supply5 = .5 will threshold the supply5
% channel and identify transitions across the .5 V level.
%
% .timeSync ; Struct passed to mdaq.readBin to synchronize the DAQ clock
%               with the neurostim clock. See mdaq.readBin
%
% BK - Oct 2023
arguments
    key % The key (ns.File and ns.ContinuousParm tuple)
    parms (1,1) struct  % The preprocessing parameters
end

%% Fetch the file to read (ns.C has already checked that it exists)
filename = fullfile(folder(ns.Experiment &key),fetch1(ns.File &key,'filename'));


%% Read using a script
tic
load(strrep(filename,'.bin','.mat'),'c');
if isfield(parms,'timeSync')
    ts = parms.timeSync;
else
    ts = struct([]); % Dont use event based timesync
end
T=  readBin(c.mdaq,file=filename,timeSync = ts);
fprintf('Read file in %d seconds.\n ',round(toc))
trialStartTime = get(ns.Experiment & key,'cic','prm','trial','what',["data" "clocktime"],'trial',1);
daqStartTime = get(ns.Experiment & key,'mdaq','prm','startDaq','what',["data" "clocktime"],'trial',1);
out = T.nsTime <= seconds(trialStartTime.clocktime(trialStartTime.data==1)/1000);
if any(out)
    fprintf('DaqStart = %s , but first trial started %d seconds later.\n',daqStartTime.data,round((trialStartTime.clocktime(trialStartTime.data==1)-daqStartTime.clocktime)/1000))
    fprintf('Removing %d samples from the mdaq that occurred before first trial start. \n',sum(out));
    T(out,:) =[];
end 

if exists(ns.Plugin & key & 'plugin_name="scanbox"')
    % Hack for a special case where mdaq is used to record scanbox
    % triggers.
    info = sbx.readInfoFile(ns.Experiment & key);    
    laserLowToHigh = find([T.laserOnDig(1);diff(T.laserOnDig)>0]);
    laserHighToLow =  find([false;diff(T.laserOnDig)<0]);
    % A ttl pulse is sent for each plane in multiplane recordings. Trim
    % down to 1 per stack. (All ttls will be stored as events below)
    laserLowToHigh = laserLowToHigh(1:info.nrPlanes:end);
    laserHighToLow = laserHighToLow(1:info.nrPlanes:end);
    nrTTL = numel(laserLowToHigh);
    if nrTTL >  info.nrFrames +1 
        % 1 extra is expected and handled in sbx.read. Try to remove
        % extraneous TTL here.
        extraTTL = nrTTL-(info.nrFrames+1);
        lastExtraTTLIx = laserHighToLow(extraTTL);
        timeToRemove = T.nsTime(lastExtraTTLIx)-T.nsTime(1);
        if timeToRemove < seconds(4)
            T(1:lastExtraTTLIx,:) =[];
        else
            error('The mismatch in TTL vs. Frames is too large to remove (%s)',timeToRemove)
        end
    elseif nrTTL < info.nrFrames
        error('TTL vs. Frames mismatch %d. This mdaq file cannot be aligned with SBX frames',nrTTL-info.nrFrames)
    end
end

%% Analog channel processing
missingAnalog  = setdiff(parms.analogChannel,T.Properties.VariableNames);
if ~isempty(missingAnalog)
    fprintf("Some analog channels were not recorded in %s. Skipping: \n",filename)
    missingAnalog %#ok<NOPRT>
end
[analogChannels] = intersect(parms.analogChannel,T.Properties.VariableNames);
time = seconds(T.nsTime);
signal = T{:,analogChannels};
[nrSamples,nrChannels] = size(signal);


%% Filter the analog channels
[signal,time] = ns.CFilter(signal,time,parms);

%% Create events based on analog or digital channels
% (Transitions between high  and low states).
if isfield(parms,'digitalChannel')
    missingDigital  = setdiff(parms.digitalChannel,T.Properties.VariableNames);
    if ~isempty(missingDigital)
        fprintf("Some digital channels were not recorded in %s. Skipping: \n",filename)
        missingDigital %#ok<NOPRT>
    end
    [digitalChannels] = intersect(parms.digitalChannel,T.Properties.VariableNames);
    if ~isempty(digitalChannels)
        % To use the same thresholding logic and code as the analog
        % channels, set asEvent threshold to 0.5
        for ch = string(digitalChannels)
            parms.asEvent.(ch) = 0.5;
        end
    end
end

% Create events for digital channels and for a selected set of analog
% channels
if isfield(parms,'asEvent')
    channelsAsEvent  = intersect(string(fieldnames(parms.asEvent)),T.Properties.VariableNames);
    plgKey  = fetch(ns.Plugin & key & 'plugin_name=''mdaq''');

    for ch = channelsAsEvent(:)'
        threshold = parms.asEvent.(ch);        
       [tf,channelIx]= ismember(ch,analogChannels);
       if tf
           % Analog channel; use the filtered signal
           thisChannel = signal(:,channelIx);
       else
           % digital channel; use the raw signal (0
           thisChannel = T{:,ch};
       end
        if isnan(threshold)
            % Estimate, assuming bimodal
            vLowHigh = (prctile(thisChannel, [5 95]));
            dV = diff(vLowHigh);
            threshold = vLowHigh(1) +dV/2; % mid point
            inMoat  = round(100*mean(thisChannel >vLowHigh(1)+dV/4 & thisChannel  < vLowHigh(2)-dV/4));
            fprintf('Analog channel (%s)  thresholding with dV %.2f:  %d%% in moat\n',ch,dV,inMoat)
        end
        digital = thisChannel>threshold;
        % Create an event representing the low to high transition of
        % the digital input
        name =ch + "High";
        goesHighNsTime = time([false; diff(digital)>0])*1000;
        fprintf('Creating %d new %sHigh events in mdaq plugin.\n',numel(goesHighNsTime), ch)
        [trialTime,trial] = eTime2TrialTime(c.mdaq.prms.vendor,goesHighNsTime);
        nrEvents= numel(goesHighNsTime);
        if nrEvents>0
            addNew(ns.PluginParameter,plgKey,name,ones(nrEvents,1),'Event',trialTime,trial,goesHighNsTime,true); % Replace existing
        end
        % Create an event representing the high to low transition of
        % the digital input
        name = ch + "Low";
        goesLowNsTime= time([false; diff(digital)<0])*1000;
        fprintf('Creating %d new %sLow events in mdaq plugin. \n',numel(goesLowNsTime), ch)
        [trialTime,trial] = eTime2TrialTime(c.mdaq.prms.vendor,goesLowNsTime);
        nrEvents= numel(goesLowNsTime);
        if nrEvents>0
            addNew(ns.PluginParameter,plgKey,name,ones(nrEvents,1),'Event',trialTime,trial,goesLowNsTime,true);% Replace existing
        end
    end
end


% Information per channel
channelInfo =struct('name',analogChannels,'nr',num2cell(1:nrChannels));
recordingInfo = struct('vendor',c.mdaq.vendor);
% MDaq has regular sampling so reduce time representation
time = [1000*time(1) 1000*time(end) nrSamples];
% Reduce storage (ns.C.align converts back to double
signal  = single(signal);

end



