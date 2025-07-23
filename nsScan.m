function [tSubject,tSession,tExperiment,isLocked] = nsScan(pv)
% This function scans filenames in (and around) the folder corresponding
% to a certain date under a given root folder and creates meta data tables
% for the subjects, sessions, and experiments found in those folders.
% Meta data are read from the file names and from JSON files.
%
% This works because Neurostim always generates files in the following
% format:
% root\YYYY\MM\dd\subject.paradigm.startTime.mat
%
% JSON files with ***additional*** meta data can be stored at three levels:
% Meta data for an experiment:
%   root\YYYY\MM\dd\subject.paradigm.startTime.json
% Meta data for a session
%   root\YYYY\MM\dd\subject.json
% Meta data for a subject
%   root\subjects.json
%
% The required metadata fields are specified in
%   root\experiment_definition.json
%   root\session_definition.json
%   root\subject_definition.json
%
% OR, if a metaDefinitionTag is specified in
% root\experiment_definition_tag.json
%   root\session_definition_tag.json
%   root\subject_definition_tag.json
% This allows different projects to use different meta data definitions in
% the same data tree.
%
% INPUT
% 'root' - Top level folder (e.g. Z:\). Folders below this root folder should
%           correspond to the years. This defaults to the value in
%           getenv('NS_ROOT')
%
% date - Date to scan.  [today]
% schedule - 'y' - Scan the year in which the date falls.
%          - 'm' - Scan the month
%          - ['d'] - Scan the specified day only.
%           - ["2022" "2023"] - scans these two years
% readJson - Read the json files containing meta data and add to the tables
%               [true]
% metaDefinitionTag - Determines which meta definition files to use (e.g. for
% metaDefinitionTag=  'rvs', the scan looks for root\experiment_definition_rvs.json
% for experiment meta data definitions.  Defaults to '' for global
% defaults in experiment_definition.json
%
% paradigm = Cell array or char of paradigms to include. Leave empty to include
%               all. Paradigms are matched case-insensitively. [{}]
%                A better way to define which paradigms to include is to
%                use the ns.Paradigm table.
% subject = Cell array or char of subjects to include. Leave empty to include all.
%               Subjects are matched case insensitivelly. [{}].
% folderFun -  Function, provided by the user, to handle special folders.
%           The nsScan function will call this function with a list of
%           folders. The user is responsible for returning a table from
%           the relevant folders. See nsScanDicomFolder for an example. Not
%           commonly needed. [].
%
% excludeSubject - Cell array of subjects to exclude. Defaults to  {'0'}
% which is the default test subject.
%
% The nsAddToDataJoint function can add the output of nsScan to a DataJoint
% database. To call that automatically, use the following two options:
%
% addToDataJoint - Set to true to add the new
%                   Subjects/Sessions/Experiments/Files to the DataJoint
%                   database. [false]
% newOnly       - Set to true to process only experiments that are not
%               already in Datajoint
% jsonOnly  - Set to true to process only the JSON meta data and update the
%               Datajoint database with those entries. (This option will
%               overrule readFileContents to false)
% minNrTrials  - include only experiments with at least this number of
% trials (also sets readFileContents =true).
% readFileContents - Set to true to read file contents when adding to the
%                   DataJoint database [false]
% cicOnly   - Set tot true to only add CIC information (and not other
%               plugins)
% safeMode  -   Set to false to remove tuples from the Datajoint database
%               without asking for confirmation (i.e., when updating information).
% OUTPUT
% tSubject      - Subject level meta data table
% tSession      - Session level meta data table
% tExperiment   - Experiment level meta data table
% isLocked      - Struct with one logical per table to indicate whether
%                   additional meta data fields are allowed or not. This is
%                   set by the Lock property in the _definition.json.
%
%
% EXAMPLES
% Scan for new files for a list of paradigms  according to a schedule
% nsScan('date',date,'schedule',schedule,'paradigm',paradigms, ...
%                    'readFileContents',true,'addToDataJoint',true,'safeMode',false, ...
%                    'newOnly',newOnly,'minNrTrials',minNrTrials);
%
% Update JSON data for existing data  (After adding manual info to the JSON
% in nsMeta)
% nsScan('date',date,'schedule',schedule,'paradigm',paradigms, ...
%                    'addToDataJoint',true,'safeMode',false, ...
%                    'jsonOnly',true);

% BK - Jan 2023
arguments
    pv.root (1,1) string = getenv('NS_ROOT')
    pv.date (1,1) datetime = datetime('now')
    pv.schedule (1,:) string = "d"
    pv.readJson (1,1) logical = true
    pv.paradigm (1,:) string = ""
    pv.subject (1,:) string =""
    pv.excludeSubject (1,:) string = "0"
    pv.folderFun = []
    pv.addToDatajoint (1,1) logical = false
    pv.newOnly (1,1) logical = false
    pv.jsonOnly (1,1) logical = false
    pv.readFileContents (1,1) logical = false
    pv.populateFile (1,1) logical = true
    pv.cicOnly (1,1) logical = false
    pv.safeMode (1,1) logical = true
    pv.metaDefinitionTag (1,1) string = ""
    pv.minNrTrials (1,:) double = 1
    pv.verbose (1,1) logical = true
    pv.fileType (1,1) string = "*.mat"
    pv.analyze  (1,1) logical = true % Show files only when analyze= true
end
assert(exist(pv.root,"dir"),"The root folder (%s) does not exist.",pv.root);

tExperiment = table;
tSession = table;
tSubject = table;
isLocked= struct('experiment',false,'session',false,'subject',false);
if pv.metaDefinitionTag~="" &&  ~startsWith(pv.metaDefinitionTag,"_")
    pv.metaDefinitionTag = "_" + pv.metaDefinitionTag;
end

%% A. Find relevant files
switch (pv.schedule)
    case "y"
        % Allow more than one folder below year
        srcFolder = fullfile(pv.root,char(datetime(pv.date,'Format','yyyy')),'**');
    case "m"
        % Exactly one folder below month
        srcFolder = fullfile(pv.root,char(datetime(pv.date,'Format','yyyy/MM')),'**');
    case "d"
        % Inside the day folder
        srcFolder = fullfile(pv.root,char(datetime(pv.date,'Format','yyyy/MM/dd')));
    case "a"
        % All year folders below the root
        % Find the year folders
        yearFolders = dir(pv.root);
        stay = regexp({yearFolders.name},'^\d{4}$','once');
        yearFolders(cellfun(@isempty,stay)) = [];
        tSubject =table;
        tSession = table;
        tExperiment =table;
        % Then a recursive call to this function for each year.
        args = namedargs2cell(pv);
        for yr = 1:numel(yearFolders)
            year = yearFolders(yr).name;
            [thisTSubject,thisTSession,thisTExperiment,isLocked] = nsScan(args{:},'schedule','y','date',[year '-01-01']);
            tSubject = [tSubject;thisTSubject]; %#ok<AGROW>
            tSession = [tSession;thisTSession];%#ok<AGROW>
            tExperiment = [tExperiment;thisTExperiment];%#ok<AGROW>
        end
        return    
    otherwise
        % Interpreting the pv.schedule as a vector of years.
        tSubject =table;
        tSession = table;
        tExperiment =table;
        % Then a recursive call to this function for each year.     
        args = namedargs2cell(pv);
        for yr = pv.schedule            
            [thisTSubject,thisTSession,thisTExperiment,isLocked] = nsScan(args{:},'schedule','y','date',yr + "-01-01");
            tSubject = [tSubject;thisTSubject]; %#ok<AGROW>
            tSession = [tSession;thisTSession];%#ok<AGROW>
            tExperiment = [tExperiment;thisTExperiment];%#ok<AGROW>
        end
        return
end

if pv.verbose
    fprintf('Scanning %s ...\n',srcFolder)
end

% Find the files matching the wildcard
dirInfo= dir(fullfile(srcFolder,pv.fileType));

if isempty(dirInfo)
    if pv.verbose
        fprintf('No %s files found in this folder: (%s)\n Is the root folder (%s ) correct? \n',pv.fileType,srcFolder,pv.root);
    end
    return;
end
[~,~,fileExtension]= fileparts(dirInfo(1).name);
% Select files matching Neurostim format filename
fullName = fullfile({dirInfo.folder}',{dirInfo.name}');
if strcmpi(filesep','\')
    fs = '\\';
else
    fs = filesep;
end
% The first part of the pattern should be thre root. But on some HPC
% systems searching in one location returns files from a different location
% Hence we allow the initial match to anything \w\d/.
pattern = ['[\w\d' fs ']*(?<session_date>\d{4,4}' fs '\d{2,2}' fs '\d{2,2})' fs '(?<subject>\w{1,10})\.(?<paradigm>\w+)\.(?<starttime>\d{6,6})\' fileExtension];
meta = regexp(fullName,pattern,'names');
% Prune those file that did not match
out = cellfun(@isempty,meta);
meta =[ meta{~out}];
dirInfo(out) =[];
fullName(out) = [];
if isempty(fullName)
    if pv.verbose
        fprintf('No Neurostim data files found in (%s)\n. Is the root folder (%s ) correct? \n',srcFolder,pv.root);
    end
    return;
end

%% Extract meta data from the file name
% and match the formats stored in the DataJoint database
[~,file,ext] = fileparts(fullName);  % Store only the filename. Path can be reconstructed from date and root.
file = strcat(file,ext);
if ~iscell(file);file={file};end
[meta.file] = deal(file{:});
[meta.bytes] = deal(dirInfo.bytes);
[meta.provenance] =deal("");
session_dates = deal(strrep({meta.session_date},filesep,'-'));  % Match ISO format of DJ
[meta.session_date] = deal(session_dates{:});
starttimes = cellfun(@(x)([x(1:2) ':' x(3:4) ':' x(5:6)]),{meta.starttime},'uni',false);
[meta.starttime] = deal(starttimes{:}); % Match the  HH:MM:SS format of DJ

if ~isempty(pv.folderFun)
    % Special handling of sub folders. Not usually needed, but see nsScanDicomFolder
    % for an example of how to include DICOM MRI data folders.
    % Note that all files in one of these folders (i.e. those returned by
    % folderFun) will be added to the ns.Files table.
    folders= dir(fullfile(srcFolder,'*')); %
    folders(~[folders.isdir])=[];  % Remove non folders
    folderFullName = strcat({folders.folder}',filesep,{folders.name}');
    metaFromFolderFun = pv.folderFun(folderFullName);
    if ~isempty(metaFromFolderFun)
        meta = catstruct(2,meta,metaFromFolderFun);
    end
end

% include specific subjects (or all if pv.subject== "")
stay = true(size(meta));
if pv.subject ~=""
    stay = stay &  ismember(upper({meta.subject}),upper(pv.subject));
end
% exclude a specific subject
stay = stay &  ~ismember({meta.subject},pv.excludeSubject);
% include based on paradigm
if pv.paradigm==""
    % Select on the basis of the ns.Paradigm table
    [pv.paradigm,pv.minNrTrials,from,to]= fetchn(ns.Paradigm,'name','mintrials','from','to');
    pv.paradigm= upper(string(pv.paradigm));
else
    % Use the specified list of paradigms to select
    pv.paradigm = upper(pv.paradigm);
    if isscalar(pv.minNrTrials)
        pv.minNrTrials = pv.minNrTrials*ones(size(pv.paradigm));
    else
        assert(numel(pv.minNrTrials)==numel(pv.paradigm),"minNrTrials should match the paradigms (or be a scalar)");
    end
    from = cell(size(pv.paradigm));
    to = cell(size(pv.paradigm));
end
if ~isempty(pv.paradigm)
    % Remove non-matching based on paradigm and (for ns.Paradigm based selection) from/to
    pdmMatch = false(size(stay));
    for pdm=1:numel(pv.paradigm)
        thisMatch = ~cellfun(@isempty, regexpi({meta.paradigm},pv.paradigm{pdm}));
        if ~isempty(from{pdm})
            thisMatch =thisMatch &  cellfun(@(x) x>=from(pdm),{meta.session_date});
        end
        if ~isempty(to{pdm})
            thisMatch =thisMatch &  cellfun(@(x) x<=to(pdm),{meta.session_date});
        end
        pdmMatch = pdmMatch | thisMatch;
    end
    stay = stay & pdmMatch;
end


meta = meta(stay);
fullName=fullName(stay);
nrExperiments = numel(meta);

if pv.newOnly || pv.jsonOnly
    if pv.jsonOnly
        % Update json based meta information for existing files only
        stay = false(1,nrExperiments);
        for i=1:nrExperiments
            stay(i) = exists(ns.Experiment & meta(i));
        end
    elseif pv.newOnly
        % Ignore experiments that are already in datajoint
        stay = true(1,nrExperiments);
        for i=1:numel(meta)
            stay(i) = ~exists(ns.Experiment & meta(i));
        end
    end
    meta = meta(stay);
    fullName=fullName(stay);
    nrExperiments = numel(meta);
end

if  nrExperiments> 0 && (pv.readFileContents || any(pv.minNrTrials >1))  && ~pv.jsonOnly
    % Read meta data from the content of the Neurostim files
    % tmp has the meta data, c the CIC objects which are passed to
    % nsAddToDataJoint below
    c = cell(1,nrExperiments);
    tmpMeta = cell(1,nrExperiments);
    out =false(1,nrExperiments);

    for i=1:nrExperiments
        try
            [tmp,c{i}] = ns.Experiment.readCicContents(meta(i),'root',pv.root);
            tmpMeta{i} = mergestruct(meta(i),tmp); % Merge to keep json/provenance meta.
        catch
            if pv.verbose
                fprintf(2,'Failed to read %s. Skipping\n',meta(i).file);
            end
            lasterr
            out(i)=true;
        end
        % Find which minium number of trials to apply (paradigm dependent)
        thisMinNrTrials = pv.minNrTrials(upper(c{i}.paradigm) == pv.paradigm);
        % Remove if too few trials
        tooFew = c{i}.nrTrialsCompleted < thisMinNrTrials;
        if tooFew && pv.verbose
            fprintf('Skipping %s ( %d trials)\n',meta(i).file,c{i}.nrTrialsCompleted);
            out(i) =true;
        end
    end
    meta  =[tmpMeta{~out}];
    c = [c{~out}];
    fullName(out) =[];
    nrExperiments = numel(meta);
else
    c= [];
end

if nrExperiments ==0
    if pv.verbose
        fprintf('No relevant files found in %s\n',srcFolder);
    end
    return;
end

%% A Experiments
tExperiment = struct2table(meta,'AsArray',true);
tExperiment = convertvars(tExperiment,@(x)(ischar(x) | iscellstr(x)),'string'); %#ok<ISCLSTR>
% Store folder root and date in the properties so that we
% can easily reconstruct filenames later
tExperiment =addprop(tExperiment,{'root'},{'table'}) ;
tExperiment.Properties.CustomProperties.root = pv.root;
if pv.readFileContents
    tExperiment = addvars(tExperiment,c','NewVariableNames','cic');
end
% Assign a unique ID to each row, so that we can match up with the original values after the user resorts.
tExperiment = addvars(tExperiment,strcat('#',string((1:nrExperiments)')),'NewVariableNames','id','before',1);
tExperiment = movevars(tExperiment,{'session_date','file','starttime','bytes','subject','paradigm'},'After',1);
metaDescriptions = ["","Date of the Session (yyyy-mm-dd: ISO 8601)","Neurostim file name","Start time of the experiment (hh:mm:ss)","bytes in the file","Unique Subject ID", "Paradigm name","Provenance information"];
if (pv.readFileContents || any(pv.minNrTrials >1))
    %Already has information from the file contents (minNrTrial>0
    %orreadFileCOntents =true)
    metaDescriptions = [metaDescriptions "Number of stimuli used" "Number of blocks" "Number of conditions" "Number of trials completed" "Matlab Version" "PTB Version" "NS version" "Run" "Sequence"];
end
if pv.readFileContents
    metaDescriptions = [metaDescriptions "CIC"];
end
tExperiment.Properties.VariableDescriptions = metaDescriptions;
if pv.readJson || pv.jsonOnly
    %% Read Meta information definitions and data from JSON files.
    allMetaFromJson = [];
    % Read the current meta definition for experiments
    definitionFile = fullfile(pv.root,"experiment_definition" +  pv.metaDefinitionTag + ".json");
    if exist(definitionFile,'file')
        json = readJson(definitionFile);
        isLocked.experiment = json.Locked =="1";
        experimentMetaFields = fieldnames(json.Fields);
        experimentMetaDescription = string(struct2cell(structfun(@(x)(string(x.Description)),json.Fields,'uni',false)));
        if any(ismember({'starttime'},experimentMetaFields))
            error('The experiment meta definition file %s should only define new meta properties, not ''starttime''',definitionFile);
        end
        nrExperimentMetaFields= numel(experimentMetaFields);
    elseif pv.metaDefinitionTag==""
        % Use no meta definition
        isLocked.experiment = false;
        nrExperimentMetaFields = 0;
    else
        error('Meta definition file %s not found',definitionFile)
    end

    emptyInit = repmat("",[height(tExperiment) 1]);
    if nrExperimentMetaFields>0
        emtpyInitFromDef = repmat({emptyInit},[1 nrExperimentMetaFields]);
        tExperiment = addvars(tExperiment,emtpyInitFromDef{:},'NewVariableNames',experimentMetaFields);
        tExperiment.Properties.VariableDescriptions(end-nrExperimentMetaFields+1:end) = experimentMetaDescription';
    end


    % Read associated JSON files (if they exist)
    for i=1:nrExperiments
        jsonFile= regexprep(fullName{i},['(\' fileExtension '$)'],'.json'); % Swap extension
        if exist(jsonFile,'file')
            thisJson = readJson(jsonFile);
        else
            thisJson = struct;
        end
        metaFieldsFromJson = fieldnames(thisJson);
        for j=1:numel(metaFieldsFromJson)
            tExperiment{i,metaFieldsFromJson{j}} = thisJson.(metaFieldsFromJson{j});
        end
        % If a definition file is used the default meta is "" but if this
        % is meta data read from a json file without a definition, then
        % every meta without a json will be initialized as []. Here we
        % replace those with  ""
        % Collect all meta
        if isempty(allMetaFromJson)
            allMetaFromJson = string(metaFieldsFromJson);
        else
            allMetaFromJson = unique([allMetaFromJson;string(metaFieldsFromJson)]);
        end
    end
    % Make sure undefined meta are ""
    for i=1:nrExperiments
        for m = allMetaFromJson'
            if isempty(tExperiment(i,m))
                tExperiment{i,m} = "";
            end
        end
    end

end
% Sort to get a consistent order
tExperiment =  sortrows(tExperiment,{'starttime','subject','paradigm'},'ascend');
if ismember("analyze",tExperiment.Properties.VariableNames) && pv.analyze
    % Remove experiments that are marked analyze = 0 in the json meta data file
    tExperiment(tExperiment.analyze=="0",:)=[];
    if height(tExperiment)==0
        fprintf('No relevant Neurostim experiment files are marked for analysis\n');
        return;
    end
end
clear meta fullName %  These are sorted differently; prevent accidental use below.

nrExperiments = height(tExperiment);    
if pv.verbose
     fprintf('Found %d matching Neurostim files: \n',nrExperiments)
     tExperiment %#ok<NOPRT>
end


%% B. Session Meta Data
% Retrieve definition and initialize table
tSession = unique(tExperiment(:,{'session_date','subject'}));
tSession.Properties.VariableDescriptions = ["Date of the Session (yyyy-mm-dd: ISO 8601)","Subject ID"];
nrSessions = height(tSession);
if pv.readJson || pv.jsonOnly
    definitionFile = fullfile(pv.root,"session_definition"+ pv.metaDefinitionTag +".json");
    if exist(definitionFile,'file')
        json = readJson(definitionFile);
        isLocked.session = json.Locked =="1"; % The _definition file states that there should be no other meta fields.
        sessionMetaFields = fieldnames(json.Fields);
        sessionMetaDescription = string(struct2cell(structfun(@(x)(string(x.Description)),json.Fields,'uni',false)));
        if any(ismember({'subject','session_date'},sessionMetaFields))
            error('The session meta definition file %s should only define new meta properties, not ''session_date''',definitionFile);
        end
        nrSessionMetaFields= numel(sessionMetaFields);
    elseif pv.metaDefinitionTag==""
        nrSessionMetaFields= 0;
    else
        error('Meta definition file %s not found',definitionFile)
    end
    emptyInit = repmat("",[height(tSession) 1]);
    if nrSessionMetaFields>0
        emtpyInitFromDef = repmat({emptyInit},[1 nrSessionMetaFields]);
        tSession = addvars(tSession,emtpyInitFromDef{:},'NewVariableNames',sessionMetaFields);
        tSession.Properties.VariableDescriptions(end-nrSessionMetaFields+1:end) = sessionMetaDescription;
    end

    % See if there are any Session JSON files and put information in the table
    for i=1:nrSessions
        if nrSessions==1
            sessionJsonFile = fullfile(pv.root,strrep(tSession.session_date,'-',filesep),tSession.subject + ".json");
            % Check to see if there are any .txt session notes (subjectNr.txt)
            txtFile = fullfile(pv.root,strrep(tSession.session_date,'-',filesep),tSession.subject + ".txt");
        else
            sessionJsonFile = fullfile(pv.root,strrep(tSession.session_date{i},'-',filesep),[tSession.subject{i} '.json']);
            % Check to see if there are any .txt session notes (subjectNr.txt)
            txtFile = fullfile(pv.root,strrep(tSession.session_date{i},'-',filesep),[tSession.subject{i}  '.txt']);
        end

        if exist(sessionJsonFile,'file')
            thisJson = readJson(sessionJsonFile);
        else
            thisJson = struct;
        end
        if exist(txtFile,"file")
            txt =string(fileread(txtFile));
            if isfield(thisJson,'comments')
                if ~any(contains(thisJson.comments,txt))
                    thisJson.comments = thisJson.comments + " " + txt;
                end
            else
                thisJson.comments = txt;
            end
        end
        metaFieldsFromJson = fieldnames(thisJson);
        for j=1:numel(metaFieldsFromJson)
            if ~ismember(metaFieldsFromJson{j},tSession.Properties.VariableNames)
                tSession =addvars(tSession,emptyInit,'NewVariableNames',metaFieldsFromJson{j});
            end
            tSession{i,metaFieldsFromJson{j}} = thisJson.(metaFieldsFromJson{j});
        end
    end
end
% Sort and assign a unique ID to each row.
tSession = sortrows(tSession,{'session_date','subject'},'ascend');
tSession = addvars(tSession,strcat('#',string((1:nrSessions)')),'NewVariableNames','id','before',1);
tSession = movevars(tSession,{'session_date','subject'},'after',1);


%% C. Subject meta data
% Both the definition and the actual subject meta data are
% in the Root folder.
% Retrieve definition
tSubject= unique(tSession(:,'subject'));
tSubject.Properties.VariableDescriptions = "Unique Subject ID";
nrSubjects = height(tSubject);
if pv.readJson || pv.jsonOnly
    definitionFile = fullfile(pv.root,"subject_definition" +  pv.metaDefinitionTag + ".json");
    if exist(definitionFile,'file')
        json = readJson(definitionFile);
        isLocked.subject = json.Locked =="1"; % No new meta fields allowed according to _definition
        subjectMetaFields = fieldnames(json.Fields);
        subjectMetaDescription = string(struct2cell(structfun(@(x)(string(x.Description)),json.Fields,'uni',false)));
        if ismember('subject',subjectMetaFields)
            error('The subject meta definition file %s should only define new meta properties, not ''subject''',definitionFile);
        end
        nrSubjectMetaFields= numel(subjectMetaFields);
    elseif pv.metaDefinitionTag==""
        nrSubjectMetaFields= 0;
    else
        error('Meta definition file %s not found',definitionFile)
    end
    emptyInit = repmat("",[height(tSubject) 1]);
    if nrSubjectMetaFields>0
        emtpyInitFromDef = repmat({emptyInit},[1 nrSubjectMetaFields]);
        tSubject = addvars(tSubject,emtpyInitFromDef{:},'NewVariableNames',subjectMetaFields);
        tSubject.Properties.VariableDescriptions(end-nrSubjectMetaFields+1:end) = subjectMetaDescription;
    end
    metaDataFile =  fullfile(pv.root,'subject.json');
    if exist(metaDataFile,'file')
        thisJson = readJson(metaDataFile);
        if isempty(thisJson)
            fprintf(2,'The json file %s returned empty.\n',metaDataFile);
        else
            metaFieldsFromJson = fieldnames(thisJson);
            thisJson = struct2table(thisJson);
            nrMetaFields= numel(metaFieldsFromJson);
            for i=1:nrSubjects
                thisSubject = tSubject.subject(i)==thisJson.subject;
                if any(thisSubject)
                    for j=1:nrMetaFields
                        if ~ismember(metaFieldsFromJson{j},tSubject.Properties.VariableNames)
                            tSubject =addvars(tSubject,emptyInit,'NewVariableNames',metaFieldsFromJson{j});
                        end
                        tSubject{i,metaFieldsFromJson{j}} = thisJson{thisSubject,metaFieldsFromJson{j}};
                    end
                end
            end
        end
    end
end
% Sort and assign a unique ID to each row.
tSubject= sortrows(tSubject,{'subject'},'ascend');
tSubject = addvars(tSubject,strcat('#',string((1:nrSubjects)')),'NewVariableNames','id','before',1);
tSubject = movevars(tSubject,{'subject'},'after',1);

if pv.addToDatajoint
    if ismember('cic',tExperiment.Properties.VariableNames)
        cic = tExperiment.cic;
        tExperiment = removevars(tExperiment,'cic');
    else
        cic =[];
    end   
    nsAddToDataJoint(tSubject,tSession,tExperiment,'cic',cic,'safeMode',pv.safeMode, ...
        'root',pv.root,'populateFile',pv.populateFile, ...
        'newOnly',pv.newOnly  && ~pv.jsonOnly);
end
