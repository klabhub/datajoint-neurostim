function failed = addFiles(tbl,pv)
% Given a table of ns.Experiment, populate ns.File with associated EGI files.
% If this fails on any of the experimtns in the tbl, their tpls are
% returned in failed.
arguments
    tbl (1,1) ns.Experiment
    pv.folder (1,1) string = "" % Look in this folder for EGI mff files/folders
end
failed = [];
% Process only the experiments that do not have MFF files associated
% already and skip "experiments" that do only impedance checking.
tbl = tbl - (ns.File & 'extension=".mff"') & 'not paradigm ="ZCheck"';
if ~exists(tbl)
    fprintf('No EGI files need to be added.\n');
end
for key = fetch(tbl,'ORDER BY session_date')'
    try
        % Loop over the table of experiments
        if pv.folder ==""
            fldr = folder(ns.Experiment &key);
        else
            fldr = pv.folder;
        end
        [~,nsFile] = fileparts(file(ns.Experiment & key));
        pdm = fetch1(ns.Experiment &key,'paradigm');
        %% Find the core MFF file linked to this experiment
        candidateMff= dir(fullfile(fldr,[key.subject '.' pdm  '*.mff']));
        coreMff = struct([]);
        for f=1:numel(candidateMff)
            % Because the EGI mff files only follow the NS convention
            % partially, need to check inside the file (for the BREC event) to
            % determine whether an MFF file belongs to a NS file.
            % (files may be read multiple times because of this).
            % Read all events, but assume a sampling rate of 1kHz, which could
            % be incorrect. Because we only use the BREC event (and not its
            % timing) this is ok. When reading the signal (ephys.egi.read) we
            % extract the actual sampling rate from the EGI bin file.
            try
                evts = mff_importevents(fullfile(candidateMff(f).folder,candidateMff(f).name), 0, 1000); % from time =0 with 1Khz sampling rate
                brec = evts(strcmpi('BREC',{evts.code})); % neurostim sends this BREC event with the ns file name
                brec = brec(1); % Sometimes there is a duplicate (with the same filename)
                egiProducingNsFile =  regexp(brec.mffkey_FLNM,[key.subject '\.' pdm '\.\d{6,6}'],'match');
                if isempty(egiProducingNsFile)
                    % No match. Check if the file was renamed with nsMeta
                    jsonFile = fullfile(fldr,nsFile + ".json");
                    if exist(jsonFile,"file")
                        json = readJson(jsonFile);
                        originalFilename = fliplr(extractBefore(fliplr(brec.mffkey_FLNM),'\'));
                        if contains(json.provenance,originalFilename)
                            % OK: this was renamed after recording
                            egiProducingNsFile = {nsFile}; % Force match below
                            fprintf(2,"This file was renamed from %s\n",brec.mffkey_FLNM);
                        end
                    else
                        % Warn and force a nonmatch.
                        fprintf(2,"Skipping unmatched MFF with internal ns file name of %s\n",brec.mffkey_FLNM);
                        egiProducingNsFile  = {'!@#!@$$'};
                    end
                end
                if contains(nsFile,egiProducingNsFile{1})
                    % Match found - no need to continue.
                    coreMff = candidateMff(f);
                    fprintf('Matching Neurostim file %s with EGI file %s\n',nsFile,egiProducingNsFile{1})
                    break;
                end
            catch me
                me.message
                fprintf(2,"Could not read events from %s\n",fullfile(candidateMff(f).folder,candidateMff(f).name))
            end
        end
        % If no core MFF can be linked, something is wrong (corrupted folder)
        if isempty(coreMff)
            fprintf(2,'No MFF file found for %s in %s\n',nsFile,fldr);
            failed = [failed;key]; %#ok<AGROW>
        else

            %% Related Impedance check files (before and after)
            % An impedance check file should have zcheck as the paradigm name
            % but there is nothing else directly linking it to the neurostim
            % experiment. Match by time.

            % TODO: how are impedances stored? info files have calibration field
            % not impedance.

            candidateMff= dir(fullfile(fldr,[key.subject '.zcheck*.mff']));
            if isempty(candidateMff)
                zBeforeFile = '';
                zAfterFile = '';
            else
                % Naming convention - limit to files for the same subject
                match = regexp({candidateMff.name},[key.subject '\.zcheck_(?<day>\d{8,8})_(?<time>\d{6,6})\.mff'],'names');
                stay = ~cellfun(@isempty,match);
                oddMff = candidateMff(~stay);
                % Warn about mff files that do not match the convention...
                for i=1:numel(oddMff)
                    fprintf(['File ' fullfile(fldr,oddMff(i).name) ' does not match the NS conventions for MFF file naming.\n'])
                end
                % Keep only files that match the convention
                mffInfo  = [match{stay}]'; % Make struct array
                candidateMff(~stay)= [] ;
                % Copy the meta information from the file name to the mffInfo
                % (used to find the nearest zcheck file).
                for i=1:numel(mffInfo)
                    fn = fieldnames(candidateMff);
                    nrFields= numel(fn);
                    for j=1:nrFields
                        mffInfo(i).(fn{j}) = candidateMff(i).(fn{j});
                    end
                end
                [~, zcheckFileOrder]= sort({mffInfo.time});
                mffInfo = mffInfo(zcheckFileOrder);
                nrZ = numel(mffInfo);
                hasImpedanceInfo = nan(1,nrZ);
                for zFileCntr = 1:nrZ
                    tempImpedanceFile = mff_importinfon(fullfile(mffInfo(zFileCntr).folder, mffInfo(zFileCntr).name),1);
                    hasImpedanceInfo(zFileCntr) = isfield(tempImpedanceFile,'impedance') || isfield(tempImpedanceFile,'calibration');
                end
                if sum(hasImpedanceInfo)<numel(hasImpedanceInfo)
                    fprintf( '%d out of %d  zcheck files do not have impedance info. Those will be ignored.\n',sum(~hasImpedanceInfo), numel(hasImpedanceInfo));
                end
                mffInfo(~hasImpedanceInfo) = [];
                % Find one before, one after.
                zBeforeFile = ephys.egi.pickOne(mffInfo,key,'BEFORE');
                zAfterFile = ephys.egi.pickOne(mffInfo,key,'AFTER');
            end

            %% GPS related files
            gpsFileDir = dir(fullfile(fldr,[key.subject '*.gpsr']));
            %select one gpsFile
            if length(gpsFileDir) >1  %choose the latest one
                [~, keepWhichGPS] = max(cellfun(@(x) datetime(x,"InputFormat","dd-MMM-uuuu HH:mm:ss"),{gpsFileDir.date}));
                gpsFileDir = gpsFileDir(keepWhichGPS);
            end
            coordinatesFileDir = dir(fullfile(fldr,[key.subject '.coordinates.sfp']));
            if isempty(coordinatesFileDir) && ~isempty(gpsFileDir)
                fprintf('The GPS file has not been solved yet - no coordinates file\n');
            elseif ~isempty(coordinatesFileDir) && isempty(gpsFileDir)
                fprintf('The GPS file does not exist. Standard template will be used for source modeling.\n');
            end
            %check if there is a file that indicates that some electrodes were
            %excluded by visual expection
            flaggedElectrodes = dir(fullfile(fldr,[key.subject '.badElectrodes.xlsx']));
            if isempty(flaggedElectrodes)
                fprintf('No electrodes were flagged as problematic through visual inspection\n');
            end


            %% Combine files, check what is already there.
            relatedFiles  = {coreMff,zBeforeFile,zAfterFile,gpsFileDir,coordinatesFileDir,flaggedElectrodes};
            out = cellfun(@isempty,relatedFiles); %Remove files that are not found
            relatedFiles(out) =[];
            relatedFiles = catstruct(1,relatedFiles{:});
            filesInDJ= ns.File & key & struct('filename',{relatedFiles.name});
            if exists(filesInDJ)
                [~,ix]= setdiff({relatedFiles.name},{fetch(filesInDJ,'filename').filename});
                relatedFiles = relatedFiles(ix);
            end

            if ~isempty(relatedFiles)
                %% Add files inside the .mff "file" ( a zipped folder)
                mffFiles = find(endsWith({relatedFiles.name},'.mff','IgnoreCase',true));
                for f = mffFiles
                    thisMff = dir(fullfile(relatedFiles(f).folder,relatedFiles(f).name));
                    thisMff([thisMff.isdir]) = [];
                    relatedFiles = catstruct(1,relatedFiles, thisMff);
                end
                [relatedFiles.isdir] = deal(false); % MFF are tagged as dirs (they are zipped dirs), but ns.File will skip those. Tag as not dir.
                if pv.folder~=""
                    %Hack- during a session the EGI files are on the iMAC,
                    %but the datajoint setup assumes they are in the folder
                    %where the experiment file lives. 
                    % pv.folder identifies the folder where mffs live on
                    % the imac.
                    % Here we create a symbolic link from the neurostim
                    % computer to the egi computer to trick datajoint in
                    % loading the correct data.
                    % This requires developer mode on windows (for the
                    % mklink /d)
                    exptFldr = folder(ns.Experiment &key);
                    for i=1:numel(relatedFiles)
                        if endsWith(relatedFiles(i).name,'.mff')
                            lnk = fullfile(exptFldr,relatedFiles(1).name);
                            if ~exist(lnk,"dir")
                                trg = fullfile(relatedFiles(1).folder,relatedFiles(1).name);
                                cmd = sprintf("mklink /d %s %s",lnk,trg);
                                [status,msg] = system(cmd);
                                if status~=0
                                    msg
                                    error('Failed to create a symlink to the EGI computer.')                                   
                                end
                            end
                        end
                        subFolder =extractAfter(relatedFiles(i).folder,pv.folder);
                        relatedFiles(i).folder =fullfile(exptFldr,subFolder);
                    end
                end
                updateWithFiles(ns.File,key, relatedFiles);
            end
        end
    catch me
        fprintf(2,"%s \n",me.message)
        fprintf(2,"Skipping. ")
    end
end

end