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
% Neurostim file. Here we find the nearest one and check inside to confirm
% that this is indeed the matching file.
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
% This function creates the following evets in s.info.egi:
%  1. All TCP Events (BREC,EREC, BTRL, ETRL, as logged by EGI)
%  2. BTRLNS - Logs when the BTRL code was sent to EGI (logged in NS)
% 3. All Diode DIN events (As logged by EGI)
%  4. Additional information calculated in this read file:
%      mffFile: the mffFile with the EEG data for this experiment
%       clockFit:  linear fit parameters based on all TCP events
%       btrlFit: linear fit parameters based on BTRL event.
%       zBeforeFile: mff file with impedance check before this  experiment
%       zAfterFile:  mff file with impedance chech after this  experiment
%       gpsFile:        GPS file for this experiment (same day).
%       coordinatesFile: Solved coordinates file based on GPS
%       dinOffsetMean:   Mean offset in ms between DIN events and TCP events.
%       dinOffsetStd:       Standard deviation of offset between DIN and TCP events.
%
%  BK - Jul 2019

%% Fetch the file to read (ns.C has already checked that it exists)


%% Read events
BEGINTIME = 0;
EGIEVENT_SAMPLINGRATE =1000;
tmp  = mff_importinfo(mffFilename); % Used to determine start time of the mff file.
fileStartTime = datetime(tmp.recordtime(1:26),'InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
evts = mff_importevents(mffFilename, BEGINTIME, EGIEVENT_SAMPLINGRATE); % from time =0 with 1Khz sampling rate
code = {evts.code};
source= {evts.sourcedevice};

brec = evts(strcmpi('BREC',code)); % neurostim sends this BREC event
if isempty(brec)
    siberr(['No BREC event found in ' mffFile.name '. Cannot match this EGI file to Neurostim']);
else
    evts(strcmpi('BREC',code)).mffkey_TRIA= '1'; % Force it to be in TRIAL 1 (not defined)

    [~,nsFile,~] = fileparts(s.file);
    [~,mffFile,~] = fileparts(mffFilename);
    if ~contains(brec.mffkey_FLNM,nsFile)
        % This MFF file was not created by the current neurostim file.
        % Maybe renamed? An error seems appropriate...but we may need a
        % way to issue a warning instead and proceed regardless
            siberr(['The MFF file was created by a different Neurostim file: ' brec,nffkey_FLNM]);
    end
    
    % Determine the time of all events (on the EGI clock)
    egiTime = char(evts.begintime);
    egiTime = datetime(egiTime(:,1:26),'InputFormat','uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
    egiTime  = seconds(egiTime - fileStartTime);  % Time in seconds since start of EGI file.
    
    %% Each of the TCP events originates from Neurostim (and should have a counterpart in the
    % s.info.egi.eventCode event). Here we store each of the events
    % that EGI registered in the sib object by using the TRIALTIME
    % (passed by neurostim to EGI). Each event also stores the time at
    % whcih EGI registered the event (in seconds since start of file).
    
    s=align(s,s.info.recon);
    
    isTcp = strcmpi(source,'Multi-Port ECI');
    tcpEvents = evts(isTcp);
    id =  code(isTcp);
    trial = str2num(fillmissing(char(tcpEvents.mffkey_TRIA),'previous')); %#ok<ST2NM> % Trial counter is update with BTRL events and applies to all subsequent events until the next update.
    targetFillConstant = size(char(tcpEvents.mffkey_TTIM),2)-4;
    fillInSpaces = repmat(' ',[1 targetFillConstant]); %occasionally, the necessary number spaces differs
    trialTime = str2num(fillmissing(char(tcpEvents.mffkey_TTIM),'constant',['-Inf' fillInSpaces])); %#ok<ST2NM> % We place all events without neurostim trial time at -inf trial time.
    data    = egiTime(isTcp); % Seconds since start of file.
    
    [s.info,~,egiSrcNr] = checkaddsource(s.info,'egi');
    [s.info] = add(s.info,egiSrcNr,id,trialTime',trial',data'); % Add all tcp events
    
    
    s.sortInfo.trialStart = s.info.egi.BTRL.data; % This is time in egi clock
    s.sortInfo.trialStop = [s.sortInfo.trialStart(2:end) egiTime(end)];
    s.sortInfo.nsTrialStart = retrieve(s.info,'clocktime',s.info.neurostim.firstframe);     % Taking first frame because .trial could be an outlier in the first trial (eye calibration).
    s.sortInfo.offset= zeros(size(s.sortInfo.trialStart)); % not used
    
    
    % Compare the time as logged by Neurostim and the time logged by
    % EGI to detect any clock drift and/or offset.
    uTcpEvents = unique(id);
    nrUEvents = numel(uTcpEvents);
    mffTime = [];
    nsTime = [];
    nsTrial = [];
    % Note that the trialtime (TTIM) stored in the egi file is actually the
    % time when the events were sent, not when they occurred. Because we
    % queue events during the trial (to avoid framedrops), and only send
    % them in the afterTrial function, most of these events will have a
    % trialTime that is very late in the trial. Of course this means they
    % are not particularly useful for alignment. We can still use them to
    % check clock drift and we do so beow.
    for i =1:nrUEvents
        mffTime = cat(2,mffTime,1000*s.info.egi.(uTcpEvents{i}).data);
        nsTime   = cat(2,nsTime,s.info.egi.(uTcpEvents{i}).clocktime);
        nsTrial  =cat(2,nsTrial,s.info.egi.(uTcpEvents{i}).trial);
    end
    if numel(mffTime(nsTrial>1))>10
        [parms] = polyfit(mffTime(nsTrial>1),nsTime(nsTrial>1),1);
        sibwarn(['Average Clock drift is ' num2str((parms(1)-1)) ' ms and the offset is ' num2str(parms(2)) ' ms'],2);
    else
        sibwarn('Not enough events to check clock drift');
        parms = [NaN NaN];
    end
    
    % In egiReadEEGorStim, we parse the EEG samples by using the BTRL event
    % (as timed by the EGI clock), to align with neurostim , we determine
    % the corresponding eventCode event with 'BTRL' as the data.
    eventCode  = s.info.egi.eventCode.data;
    eventCodeTrialTime = s.info.egi.eventCode.time;
    eventCodeTrial = s.info.egi.eventCode.trial;
    isBtrl = strcmpi('BTRL',eventCode);
    % Add this as an event called BTRLNS
    s.info = add(s.info,egiSrcNr,'BTRLNS',eventCodeTrialTime(isBtrl)',eventCodeTrial(isBtrl)');
    
    % Here too we can check clock drift
    btrlParms=polyfit(1000*s.info.egi.BTRL.data,s.info.egi.BTRLNS.clocktime,1);
    sibwarn(['BTRL Clock drift is ' num2str((btrlParms(1)-1)) ' ms and the offset is ' num2str(btrlParms(2)) ' ms'],2);
    
    
    %% DIN1 syncs
    % Some files will have DIN1 (diode) events that are linked to neurostim
    % and TCP events. We store the din events for analysis, and we check
    % some of the alignment.
    isDin = strncmpi(code,'DIN',3);
    if any(isDin)
        % Find out who generated the flashes.
        c= s.info.cic;
        logsOnset = find(cellfun(@(x) (~isempty(x) && contains(func2str(x),'neurostim.plugins.egi.logOnset')),{c.stimuli.onsetFunction}));
        logsOffset = find(cellfun(@(x) (~isempty(x) && contains(func2str(x),'neurostim.plugins.egi.logOffset')),{c.stimuli.offsetFunction}));
        potentialLoggers = [logsOnset logsOffset];
        nrPotentialLoggers = numel(potentialLoggers);
        eventName = '';
        %make sure that we do not register a DIN event even though there
        %was none
        dinOffsetMean =NaN;
        dinOffsetStd=NaN;
        if nrPotentialLoggers ==0
            sibwarn('DIN events found in EGI file, but no stimulus has a egi.logOnset function.');            
        else
            for i=potentialLoggers       
                stimulusFields = fieldnames(c.stimuli(i));
              	if any(strcmpi(stimulusFields,'diode'))
                    if c.stimuli(i).diode.on && ~c.stimuli(i).diode.whenOff
                        eventName = [c.stimuli(i).name(1:2) 'ON'];
                    else
                    
                        if c.stimuli(i).diode.on && c.stimuli(i).diode.whenOff
                            eventName = [c.stimuli(i).name(1:2) 'OF'];
                        end
                    end
                    
                elseif any(strcmpi(stimulusFields,'diodeFlasher'))
                    if c.stimuli(i).diodeFlasher.enabled && isequal(c.stimuli(i).diodeFlasher.onColor,[1 1 1])
                        % This stimulus created a flash for the diode when it
                        % turned on. This moment was logged Therefore there should be a matching DIN event
                        % for each [stimulus.name(1:2) ON] event.
                        eventName = [c.stimuli(i).name(1:2) 'ON'];
                    else
                        if c.stimuli(i).diodeFlasher.enabled && isequal(c.stimuli(i).diodeFlasher.onColor,[0 0 0])
                            % This stimulus flashed the diode on when the stimulus went
                            % off. So there should be prefixOF events matching DIN
                            eventName = [c.stimuli(i).name(1:2) 'OF'];
                        end
                    end
                end
                
                if isempty(eventName)
                    sibwarn(['egi.log onset/offset function, but no diode on in ' c.stimuli(i).name]);
                    continue;
                else
                    isMatchedTcp = strcmpi(code,eventName);
                end
                dinEvts =evts(isDin);
                tcpEvts = evts(isMatchedTcp);
                tDin =egiTime(isDin);
                tTcp = egiTime(isMatchedTcp);
                nrDin = numel(tDin);
                nrTcp = numel(tTcp);
                if nrDin~=nrTcp
                    dt = finddelay(tDin,tTcp); % Negative means tDin has some extra events at the beginning
                    if dt<0
                        dinEvts = dinEvts(abs(dt)+1:end);
                        tDin = tDin(abs(dt)+1:end);
                    elseif dt>0
                        tcpEvts = tcpEvts(dt+1:end);
                        tTcp = tTcp(dt+1:end);
                    end
                end
                nrDin = numel(tDin);
                nrTcp = numel(tTcp);
                [dinOffsetMean,dinOffsetStd] = mstd(1000*(tDin-tTcp));% ms
                sibwarn(['Diode events are offset by ' num2str(dinOffsetMean,3) ' ms with a standard deviation of ' num2str(dinOffsetStd,3) ' ms'],1);
                
                %Add the DIN events to the infos
                id ={dinEvts.code};
                isBtrl = strcmpi(code,'BTRL');
                btrlTime = egiTime(isBtrl);
                cntr=0;
                trial = nan(1,nrDin);
                trialTime = nan(1,nrDin);
                for t=tDin(:)'
                    cntr= cntr+1;
                    thisTrial = find(t>btrlTime,1,'last');
                    if ~isempty(thisTrial)
                        trial(cntr) = thisTrial;
                    end
                    trialTime(cntr) = t-btrlTime(thisTrial);
                end
                % Note that we cannot compare these tcp events' trial time to the
                % eventCode trial time because the latter are logged after the end
                % of the trial
                s.info = add(s.info,egiSrcNr,id,trialTime,trial);
            end
        end
    else
        sibwarn('No DIN events found in EGI file. No photodiode?')   ;
        dinOffsetMean =NaN;
        dinOffsetStd=NaN;
    end
    
    %% Add supplemental info
    % Add information about the electrodes  - hardcoded for now
    s.eInfo.electrodes = 1:257;  % Default for now.
    s.eInfo.samplingFreq = 1000; % Will adjust later if reading EEG
    s.eInfo.breakout = []; % Mayeb stim could go here.
    
    
   
    
    ePhysFileFound = true;
    
end
end

function v = filenameOnly(f)
[~,f,e] = fileparts(f);
v = [f e];
end



function s  =egiReadEEGorSTIM(s,readwhatstr,varargin)
% Read EEG/STIM data from EGI MFF files.
% Data are stored in the raw EGI MFF files 
% EEG data will be stored in the iLFP analog object
% STIM data will be stored in the iANALOG analog object

% p =inputParser;
% p.addParameter('cleanupMode',{},@iscell); % Cleanup data. Array of cleanup modes. e.g. ('lfpClip','stg','iti'}
% p.addParameter('cleanupParms',{},@iscell); % passed verbatim to cleanupcdata; should be parm/value pairs
% p.parse(varargin{:});
% 
if ~any(strcmp(readwhatstr,{'LFP','ANALOG'}))
    error('readwhatstr must be either ''LFP (=EEG)'' or ''ANALOG (=STIM)''');
end
tt=tic;
sibwarn(['Loading EGI MFF ' readwhatstr ' data ... '],2,1);
cachealign = s.info.currentAlign;
s=align(s,s.info.recon);
switch upper(readwhatstr)
    case 'LFP'
        % Read EEG Data
        mffFile = s.info.egi.mffFile.data; % This has been discovered and set by egiread.m 
        mffFilePath = fullfile(fileparts(s.file),mffFile); %but we need the full path, otherwise it will be appended to pwd and therefore result in an invalid path

        if isempty(mffFile)
            error(s,['No MFF File for ' s.file ' has been found ']);
        end
        sibwarn(['Starting to read from ' mffFile ],1);
        EEG = mff_import(mffFilePath); % Use  Arno Delorme's Toolbox
        sibwarn(['Done reading from ' mffFile ],1);
     
        %EEG.data is the raw data [nrElectrodes nrSamples]
        timeInSeconds = (0:EEG.pnts-1)/EEG.srate; % Time in egi clock     
        electrodes = s.lfpToRead;         % these the user wants to keep.             
        target = 'iLFP';
        s.eInfo.samplingFreq = EEG.srate;
        s.eInfo.electrodes  = 1:EEG.nbchan;
               
    case 'ANALOG'
        % Read STIM data  
        siberr('Not implemented yet');
end
% Align to the event created in egiread.
s = align(s,s.info.egi.BTRL);
nrTrials = s.info.nrTrials;
ana = cell(1,nrTrials);
for trial = 1:nrTrials-1
 ana{trial} = EEG.data(electrodes,timeInSeconds>= s.sortInfo.trialStart(trial) & timeInSeconds <  s.sortInfo.trialStart(trial+1),:)';
end
% Put the rest in the last trial
ana{nrTrials} = EEG.data(electrodes,timeInSeconds>= s.sortInfo.trialStart(nrTrials),:)';

% All data are now aligned to BTRL with zero offset, the corresponding event in
% neurostim time is BTRLNS (see egiread). Store 
s.(target) = set(s.(target),ana,'e',electrodes,'sf',EEG.srate,'rawAlignEvent','BTRLNS','rawAlignSource','egi','rawAlignOffset',0,'replace',false);
s = align(s,cachealign);
sibwarn(['Done loading EGI ' upper(readwhatstr) ' data (' seconds2readable(toc(tt)) ')'],2,-1);
end

