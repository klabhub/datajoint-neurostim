function addFiles(tbl)
% Given a table of ns.Experiment, populate ns.File with associated EGI files.
arguments
    tbl (1,1) ns.Experiment
end
tbl = tbl - (ns.File & 'extension=".mff"');
for key = fetch(tbl)'
    fldr = folder(ns.Experiment &key);

    %% Find candidate MFF files
    egiFilesDir = dir(fullfile(fldr,'*.mff'));

    % Naming convention
    match = regexp({egiFilesDir.name},[key.subject '\.(?<pdm>\w+)_(?<day>\d{8,8})_(?<time>\d{6,6})\.mff'],'names');
   
   
    stay = ~cellfun(@isempty,match);
    oddMff = egiFilesDir(~stay);
    % Warn about mff files that do not match the convention...
    for i=1:numel(oddMff)
        fprintf(['File ' fullfile(fldr,oddMff(i).name) ' does not match the NS conventions for MFF file naming.\n'])
    end

    mffInfo  = [match{stay}]'; % Make struct array
    dirInfo = egiFilesDir(stay); % Files that match the NS convention
    for i=1:numel(mffInfo)
        fn = fieldnames(dirInfo);
        nrFields= numel(fn);
        for j=1:nrFields
            mffInfo(i).(fn{j}) = dirInfo(i).(fn{j});
        end
    end

    % Select only MFF with same  paradigm, and day
    stay = datetime(key.session_date,"InputFormat","uuuu-MM-dd")==datetime({mffInfo.day}',"InputFormat","uuuuMMdd");
    if ~any(stay)
        fprintf('No date matches in mff files.\n');
    end
    mffInfo = mffInfo(stay);

    samePdm= strcmpi(fetch1(ns.Experiment &key,'paradigm'),{mffInfo.pdm}');
    if ~any(samePdm)
        fprintf('No paradigm matches in mff files. \n');
    end
    keepMff = ephys.egi.pickOne(mffInfo(samePdm),key,'NEAREST'); % We have now identified a file for reading
    keepMff.isdir = false; % MFF are tagged as dirs (they are zipped dirs), but ns.File will skip those. Tag as not dir.
    updateWithFiles(ns.File,key,keepMff);
    
    mffFilename = fullfile(keepMff.folder,keepMff.name);

    %% Impedance check files    
    % An impedance check file should have the same subject, but paradigm
    % should be zcheck
    [~, mffs.subject, ~] = fileparts(mffFilename);
    mffs.subject = mffs.subject(1:3);
    stay = ismember({mffInfo.subject},key.subject) & ismember({mffInfo.pdm},'zcheck');       
    zMff = mffInfo(stay);
    [~, zcheckFileOrder]= sort({zMff.time});
    zMff = zMff(zcheckFileOrder);
    nrZ = numel(zMff);
    hasImpedanceInfo = nan(1,nrZ);
    for zFileCntr = 1:nrZ
        tempImpedanceFile = mff_importinfon(fullfile(zMff(zFileCntr).folder, zMff(zFileCntr).name),1);
        hasImpedanceInfo(zFileCntr) = isfield(tempImpedanceFile,'impedance');
    end
    
    if sum(hasImpedanceInfo)<numel(hasImpedanceInfo)
        fprintf( '%d out of %d  zcheck files do not have impedance info. Those will be ignored.\n',sum(hasImpedanceInfo), numel(hasImpedanceInfo));
    end
    zMff(~hasImpedanceInfo) = [];

    zBeforeFile = pickOne(zMff,key,'BEFORE');
    zAfterFile = pickOne(zMff,key,'AFTER');
   	
    %make sure that zBeforeFile and and zAfterfFile are not the same
    if strcmpi(zBeforeFile,zAfterFile)
    %TODO    [zBeforeFile, zAfterFile] = compareTwo(mffInfo,s,stay,zBeforeFile,mffFilename);
    end
    
    %% Related GPS files

    gpsFileDir = dir(fullfile(fldr,'*.gpsr'));
    %select one gpsFile
    if length(gpsFileDir) >1  %choose the latest one
        [~, keepWhichGPS] = max(cellfun(@(x) datetime(x,"InputFormat","dd-MMM-uuuu HH:mm:ss"),{gpsFileDir.date}));
        gpsFileDir = gpsFileDir(keepWhichGPS);
    end
    
    coordinatesFileDir = dir(fullfile(fldr,[key.subject '.coordinates.sfp']));
    
    if ~isempty(coordinatesFileDir) && ~isempty(gpsFileDir)
        fprintf('The GPS file has not been solved yet - no coordinates file\n');
    elseif ~isempty(coordinatesFileDir) && isempty(gpsFile)
        fprintf('The GPS file does not exist. Standard template will be used for source modeling.\n');            
    end
    
    %check if there is a file that indicates that some electrodes were
    %excluded by visual expection
    flaggedElectrodesDir = dir(fullfile(fldr,[key.subject '.badElectrodes.xlsx']));
    
    if ~isempty(flaggedElectrodesDir)
      	fprintf('No electrodes were flagged as problematic through visual inspection\n');   
    end
    

    %% Combine files, check what is already there, add the rest.

      updateWithFiles(ns.File,key,gpsFileDir);



end
end