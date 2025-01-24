function failed = addFiles(tbl)
% Given a table of ns.Experiment, populate ns.File with associated EGI files.
% If this fails on any of the experimtns in the tbl, their tpls are
% returned in failed.
arguments
    tbl (1,1) ns.Experiment
end
failed = [];
if ~ispc
    % This ***may*** not work well on a linux system. Not clear whether that
    % could be fixed in the Java code, or that this is a big vs. small
    % endian issue. 
    fprintf(2,"I found that mmf_importinfon sometimes cannot read info1.xml files (with impedance values) on Linux systems, when it works fine on Windows.\n Running on this machine may fail.")
end
% Process only the experiments that do not have MFF files associated
% already and skip "experiments" that do only impedance checking.
tbl = tbl - (ns.File & 'extension=".mff"') & 'not paradigm ="ZCheck"';
if ~exists(tbl)
    fprintf('No EGI files need to be added.\n');
end
for key = fetch(tbl)'
    % Loop over the table of experiments
    fldr = folder(ns.Experiment &key);
    [~,nsFile] = fileparts(file(ns.Experiment & key));
    pdm = fetch1(ns.Experiment &key,'paradigm');
    %% Find the core MFF file linked to this experiment
    candidateMff= dir(fullfile(fldr,[key.subject '.' pdm  '_*.mff']));
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
        brec = evts(strcmpi('BREC',{evts.code})); % neurostim sends this BREC event  
        egiProducingNsFile =  regexp(brec.mffkey_FLNM,[key.subject '.' pdm '\.\d{6,6}'],'match');  
        assert(~isempty(egiProducingNsFile),'This MFF file has an incorrectly formatted NS file in its BREC',candidateMff(f).name);
        if contains(nsFile,egiProducingNsFile{1})
            % Match found - no need to continue.
            coreMff = candidateMff(f);
            fprintf('Matching Neurostim file %s with EGI file %s\n',nsFile,egiProducingNsFile{1})
            break;
        end
        catch me
            fprintf(2,"Could not read events from %s\n",fullfile(candidateMff(f).folder,candidateMff(f).name))
        end
    end
    % If no core MFF can be linked, something is wrong (corrupted folder)
    if isempty(coreMff)
        fprintf(2,'No MFF file found for %s in %s',nsFile,fldr);
        failed = [failed;key];
    else
    
    %% Related Impedance check files (before and after)
    % An impedance check file should have zcheck as the paradigm name
    % but there is nothing else directly linking it to the neurostim
    % experiment. Match by time.

    % TODO: how are impedances stored? info files have calibration field
    % not impedance.

    candidateMff= dir(fullfile(fldr,[key.subject '.zcheck_*.mff']));   
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
        updateWithFiles(ns.File,key, relatedFiles);
    end
    end
end

end