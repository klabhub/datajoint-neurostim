function EEG = dataset(key,pv)
% Uses EEGLAB plugins (mffmatlabio and fieldtrip) to read egi MFF
% files and returns an EEG struct. This is used to add EGI data to the
% database (ephys.egi.read),
% 
% Neurostim adds the following fields
% EEG.etc.neurostim.clockParms   - polyval(clockParms,EGI.time) -> converts EGI time to neurostim time
% EEG.etc.neurostim.pluginparameter -  the events stored in the MFF file  packaged as a struct array for easy insertion into the pluginparameter table (used by ephys.egi.read_
% EEG.etc.neurostim.expt = the Experiment tuple
%
% See also ephys.eeglab.addEvents, ephys.egi.read
arguments
    key (1,1)      % Experiment key or a keySource for the ns.C table
    pv.data (1,1) string % RAW, EMPTY, or a ctag
    pv.plg (1,:)  string  = string.empty % Plugin name 
    pv.prm (1,:) string   = string.empty % Event names
    pv.itag (1,1) string = "" % The ICA to load, identified by its itag. "" means no ICA will be loaded.
end 

if ismember(upper(pv.data),["RAW" "EMPTY"])    
    % Ignore ctag to determine ns.C key
    key = fetch(ns.File & key & 'extension=".mff"','filename');
else
    key.ctag = pv.data;
    key = fetch(ns.C & key & 'filename LIKE "%.mff"','filename');    
end
assert(~isempty(key),"This experiment does not have an associated MFF file");
mffFile = fullfile(folder(ns.Experiment &key),key.filename);
mffFile= strrep(mffFile,'\','/'); % Avoid fprintf errors
assert(exist(mffFile),"MFF file %s does not exist.",mffFile); %#ok<EXIST>


switch pv.data
    case "RAW"
        % Call mff_import directly to read everything
        eegLabSave = 0 ; % Don't save in eeglab
        correctEvents = 0; %  Don't correct events with UTF chars/
        fprintf("Using eeglab to read header and data from " +  mffFile + "...\n")
        EEG = pop_mffimport(char(mffFile),{},eegLabSave,correctEvents);
        urSrate = EEG.srate;
    otherwise        
        % Adapted code from mff_import to avoid reading the signal
        % Some pieces (that we don't currently use) are missing
        fprintf("Using eeglab to read header from " +  mffFile + "...\n")

        %  Initialize an empty standard EEGLAB structure
        EEG = eeg_emptyset();
        %  Use the ft read header v1 tool to parse metadata without reading the binary signal
        % (mff_import always reads the signal)
        mffHeader = ft_read_header(mffFile,'headerformat','egi_mff_v1');
        %  Map the XML metadata metadata onto the empty EEGLAB structure
        EEG.setname    = sprintf('%s@%sT%s',key.subject,key.session_date,key.starttime);
        EEG.srate      = mffHeader.Fs;
        EEG.nbchan     = mffHeader.nChans;
        EEG.pnts       = mffHeader.nSamples;
        EEG.trials     = 1;

        % Import channel locations from the MFF coordinates XML
        [EEG.chanlocs, EEG.ref] = mff_importcoordinates(mffFile);
        if iscell(EEG.ref)
            EEG.ref = sprintf('%s ', EEG.ref{:});
        end
        EEG.urchanlocs = EEG.chanlocs;
        if ~isempty(EEG.chanlocs)
            EEG = eeg_checkchanlocs(EEG); % put fiducials in chanfinfo
        end
        EEG=pop_chanedit(EEG, 'forcelocs',[],'nosedir','+Y');
        EEG.chaninfo.filename = 'egimff';

        % Determine when recording started
        [~, begTime]= mff_importinfo(mffFile);

        % Import the event tracks
        correctEvents=0;
        EEG.event      = mff_importevents(mffFile,begTime,EEG.srate,correctEvents);
        urSrate = EEG.srate;
        if upper(pv.data)=="EMPTY" 
           % Set the data to all-zero sparse to avoid erors on channel/time
            % selection
            EEG.data =sparse(EEG.nbchan,EEG.pnts);
        else
            %Get the data from a ctag
            cRel = ns.C & key & struct('ctag',pv.data);
            assert(exists(cRel),'No ns.C data with ctag %s found for this experiment.',pv.data);
            fprintf("Retrieving preprocessed data from ns.C (ctag=%s)\n",pv.data)
            C = fetch(ns.CChannel &cRel,'channelinfo','signal');
            EEG.chanlocs = [C.channelinfo];
            EEG.data = [C.signal]';
            [EEG.nbchan,EEG.pnts] = size(EEG.data);
            EEG.srate  = round(cRel.samplingRate);
            EEG.xmin  =0;
            EEG.xmax = (EEG.pnts-1)/EEG.srate+EEG.xmin;

            % Check if there are ICA results
            if pv.itag ~=""
                    icaKey = key;
                    icaKey.itag = pv.itag;
                    w = ns.Ica.getWeights(icaKey);  % handles both session and per-exp ICA
                    if isempty(w)
                                       fprintf("ICA with itag %s not found.\n",pv.itag);
                    elseif ~isempty(w.chanlabels)
                        % Session ICA: some channels may have been excluded
                        % during the session ICA (bad in another experiment).
                        % Find which channels are included vs excluded here.
                        allLabels = {EEG.chanlocs.labels};
                        [~, icachansind] = ismember(w.chanlabels, allLabels);
                        icachansind = icachansind(icachansind > 0);
                        assert(~isempty(icachansind), ...
                            'None of the session ICA channels found in this experiment''s chanlocs.');

                        EEG.icachansind = icachansind;
                        EEG.icasphere   = w.sphere;
                        EEG.icaweights  = w.weights;
                        EEG.icawinv     = w.winverse;
                        EEG.icaact = icaact(EEG.data(icachansind,:), ...
                        EEG.icaweights * EEG.icasphere, mean(EEG.data(icachansind,:), 2));

                        % Project excluded channels into ICA space via OLS so
                        % that pop_subcomp can clean them too.
                        % W_excl = D_excl * A' * inv(A * A')
                        % where A = icaact  [nComps x nSamps]
                        exclIdx = setdiff(1:EEG.nbchan, icachansind);
                        if ~isempty(exclIdx)
                            A  = EEG.icaact;                        % [nComps x nSamps]
                            D  = EEG.data(exclIdx, :);              % [nExcl  x nSamps]
                            % OLS: W_excl = D * A' / (A * A')
                            W_excl = (D * A') / (A * A');           % [nExcl x nComps]
                            % Augment icawinv with projected rows for excluded channels
                            winv_aug = zeros(EEG.nbchan, size(EEG.icawinv, 2));
                            winv_aug(icachansind, :) = EEG.icawinv;
                            winv_aug(exclIdx,    :) = W_excl;
                            EEG.icawinv = winv_aug;
                            % Fold sphere into weights and extend to all channels so that
                            % eeg_checkset is satisfied: size(icaweights,2) == size(icasphere,2)
                            %                            == numel(icachansind).
                            % icaweights_aug * eye * data(1:nbchan,:) reproduces icaact
                            % because the excluded-channel columns are zero.
                            weights_aug = zeros(size(EEG.icaweights, 1), EEG.nbchan);
                            weights_aug(:, icachansind) = EEG.icaweights * EEG.icasphere;
                            EEG.icaweights  = weights_aug;
                            EEG.icasphere   = eye(EEG.nbchan);
                            EEG.icachansind = 1:EEG.nbchan;
                        end
                    else
                        % Per-experiment ICA: indices stored directly
                        EEG.icachansind = w.channels;
                        EEG.icasphere   = w.sphere;
                        EEG.icaweights  = w.weights;
                        EEG.icawinv     = w.winverse;
                        EEG.icaact = icaact(EEG.data, EEG.icaweights * EEG.icasphere, mean(EEG.data, 2));
                    end                
            end
            EEG= eeg_checkset(EEG);
        end        
 end


%% Process BREC event to ensure the neurostim and mff file correspond to the same experiment
nrEvts = numel(EEG.event);
eventCode = {EEG.event.code};
brec = EEG.event(strcmpi('BREC',eventCode)); % neurostim sends this BREC event
assert(~isempty(brec),"No BREC event found in " +  mffFile + ". Cannot match this EGI file to Neurostim");
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

%% Preprocess the BTRL events to get trial and time.
isBeginTrial =strcmpi(eventCode,'BTRL');
trial = nan(nrEvts,1);
trial(isBeginTrial) = cellfun(@str2num,{EEG.event(isBeginTrial).mffkey_TRIA});
trial(1) =1; % Group events before trial 1 with trial 1
trial = fillmissing(trial,"previous");
trial =num2cell(trial);
[EEG.event.trial]=deal(trial{:});
% Determine the time of all events (on the EGI clock)
% Use the .latency field of the MFF.event, not the begintime (which can be
% offset by a few hundred ms). .latency is in samples.
eventEgiTime = ([EEG.event.latency]-1)/urSrate; % This is seconds since the start of the data in EGI
eventEgiTime  =num2cell(eventEgiTime );
[EEG.event.egitime] = deal(eventEgiTime {:});


%% Read the properties of cic and the egi plugin for this experiment
% Synchronize clocks.
% Using NTPSync results in pretty much perfectly aligned clocks (no
% drift). But we check anyway using the Begin Trial (bTRL) events.
prms  = get(ns.Experiment &key,{'cic','egi'});
trialStartTimeNeurostim  = prms.cic.trial.clocktime(2:end);%
trialStartTimeEgi = [EEG.event(strcmpi(eventCode,'BTRL')).egitime];
% The number of trials should match
assert(numel(trialStartTimeEgi)==numel(trialStartTimeNeurostim),'Number of trials mismatched in EGI and NS');
% Determine clock drift, and the offset between the first trial start event
% in neurostim and in EGI.
if numel(trialStartTimeNeurostim)>1
    EEG.etc.neurostim.clockParms = polyfit(trialStartTimeEgi,trialStartTimeNeurostim,1);
    fprintf(['Average Clock drift is ' num2str((EEG.etc.neurostim.clockParms(1)-1000)) ' ms/s and the offset is ' num2str(EEG.etc.neurostim.clockParms(2)) ' ms \n' ]);
else
    fprintf('Single trial; assuming zero clock drift.\n');
    % Assume that the trialStartTimeEgi and trialStartTimeNeurostim refer
    % to the same time and that there is no clock drift. The 1000 slope
    % takes the s from egi to ms used here.
    EEG.etc.neurostim.clockParms = [1000 trialStartTimeNeurostim]; % Assming zero drift, zero offset
end


%% Package the events in a tpl that can be inserted into the DJ PluginParameter table
isBoundary = strcmpi({EEG.event.type},'boundary');
nsEvents = EEG.event(~isBoundary);
eventNsTime = polyval(EEG.etc.neurostim.clockParms,[nsEvents.egitime]);
eventTrial = [nsEvents.trial];
eventTrialTime = eventNsTime' - trialStartTimeNeurostim(eventTrial);
eventCode = string({nsEvents.code});
uNames= unique(eventCode);
nrNames= numel(uNames);

prmTpl  = struct('property_name','','property_time',[],'property_nstime',[],'property_trial',[],'property_value',[],'property_type','Event');
prmTpl = repmat(prmTpl,[nrNames 1]);
nmCntr =0;
for nm = uNames
    nmCntr = nmCntr+1;
    prmTpl(nmCntr).property_name = char(nm);
    stay= eventCode == nm;
    prmTpl(nmCntr).property_time = eventTrialTime(stay);
    prmTpl(nmCntr).property_nstime=eventNsTime(stay);
    prmTpl(nmCntr).property_trial= eventTrial(stay);
    prmTpl(nmCntr).property_value = nsEvents(stay);
end
EEG.etc.neurostim.pluginparameter = prmTpl;
[EEG.filepath, name, ext] = fileparts(char(mffFile));
EEG.filename = [name ext];
EEG.etc.neurostim.expt = key;

if EEG.srate ~=urSrate
    for iEvent=1:length(EEG.event)
        EEG.event(iEvent).latency = round(EEG.event(iEvent).latency*(EEG.srate/urSrate));
    end    
end
%%
% Check consistency
EEG = eeg_checkset(EEG);

if ~isempty(pv.plg)
    EEG= ephys.egi.eeglabAddEvents(EEG,pv.plg,pv.prm);
end