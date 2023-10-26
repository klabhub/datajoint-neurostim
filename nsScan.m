function [tSubject,tSession,tExperiment,isLocked] = nsScan(varargin)
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
%
% readJson - Read the json files containing meta data and add to the tables
%               [true]
% metaDefinitionTag - Determines which meta definition files to use (e.g. for
% metaDefinitionTag=  'rvs', the scan looks for root\experiment_definition_rvs.json
% for experiment meta data definitions.  Defaults to '' for global
% defaults in experiment_definition.json
%
% paradigm = Cell array or char of paradigms to include. Leave empty to include
%               all. Paradigms are matched case-insensitively. [{}]
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
% BK - Jan 2023

p = inputParser;
p.addParameter('root', getenv('NS_ROOT'));
p.addParameter('date',datetime('now'));
p.addParameter('schedule','d');
p.addParameter('readJson',true);
p.addParameter('paradigm',{});
p.addParameter('subject',{});
p.addParameter('excludeSubject',{'0'});
p.addParameter('folderFun','');
p.addParameter('addToDataJoint',false,@islogical)
p.addParameter('newOnly',false,@islogical);
p.addParameter('readFileContents',false,@islogical)
p.addParameter('populateFile',true,@islogical)
p.addParameter('populateMovie',true,@islogical)
p.addParameter('cicOnly',false,@islogical);
p.addParameter('safeMode',true,@islogical);
p.addParameter('metaDefinitionTag','',@ischar);
p.parse(varargin{:});

tExperiment = [];
tSession = [];
tSubject = [];

if ~isempty(p.Results.metaDefinitionTag) && p.Results.metaDefinitionTag(1)~='_'
    metaDefinitionTag = ['_' p.Results.metaDefinitionTag];
else
    metaDefinitionTag  = ''; % Use global default definition 
end

%% A. Find relevant files
switch (p.Results.schedule)
    case 'y'
        % Allow more than one folder below year
        srcFolder = fullfile(p.Results.root,char(datetime(p.Results.date,'Format','yyyy')),'**');
    case 'm'
        % Exactly one folder below month
        srcFolder = fullfile(p.Results.root,char(datetime(p.Results.date,'Format','yyyy/MM')),'**');
    case 'd'
        % Inside the day folder
        srcFolder = fullfile(p.Results.root,char(datetime(p.Results.date,'Format','yyyy/MM/dd')));
    otherwise
        error('Unknown schedule %s',p.Results.schedule)
end

fprintf('Scanning %s for %s ...\n',srcFolder,strjoin(p.Results.paradigm,'/'))

% Find the files matching the wildcard
dirInfo= dir(fullfile(srcFolder,'*.mat'));

if isempty(dirInfo)
    fprintf('No .mat files found in this folder: (%s)\n Is the root folder (%s ) correct? \n',srcFolder,p.Results.root);
    return;
end

% Select files matching Neurostim format filename
fullName = fullfile({dirInfo.folder}',{dirInfo.name}');
if strcmpi(filesep','\')
    fs = '\\';
else
    fs = filesep;
end
if ispc
    root = strrep(p.Results.root,'/','\');
else 
    root = p.Results.root;
end
pattern = strcat(strrep(root,'\',fs),[fs '?'], ['(?<session_date>\d{4,4}' fs '\d{2,2}' fs '\d{2,2})' fs '(?<subject>\w{1,10})\.(?<paradigm>\w+)\.(?<starttime>\d{6,6})\.mat$']);
meta = regexp(fullName,pattern,'names');
% Prune those file that did not match
out = cellfun(@isempty,meta);
meta =[ meta{~out}];
dirInfo(out) =[];
fullName(out) = [];
if isempty(fullName)
    fprintf('No Neurostim data files found in (%s)\n. Is the root folder (%s ) correct? \n',srcFolder,p.Results.root);
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

if ~isempty(p.Results.folderFun)
    % Special handling of sub folders. Not usually needed, but see nsScanDicomFolder
    % for an example of how to include DICOM MRI data folders.
    % Note that all files in one of these folders (i.e. those returned by
    % folderFun) will be added to the ns.Files table.
    folders= dir(fullfile(srcFolder,'*')); %
    folders(~[folders.isdir])=[];  % Remove non folders
    folderFullName = strcat({folders.folder}',filesep,{folders.name}');
    metaFromFolderFun = p.Results.folderFun(folderFullName);
    if ~isempty(metaFromFolderFun)
        meta = catstruct(2,meta,metaFromFolderFun);
    end
end

% Select on subject/paradigm if not empty
stay = true(size(meta));
if ~isempty(p.Results.subject)
    if ischar(p.Results.subject)
        subject = {p.Results.subject};
    else
        subject =p.Results.subject;
    end
    stay = stay &  ismember(upper({meta.subject}),upper(subject));
end

stay = stay &  ~ismember({meta.subject},p.Results.excludeSubject);

if ~isempty(p.Results.paradigm)
    if ischar(p.Results.paradigm)
        paradigm = {p.Results.paradigm};
    else
        paradigm =p.Results.paradigm;
    end
    stay = stay & ismember(upper({meta.paradigm}),upper(paradigm));
end

meta = meta(stay);
fullName=fullName(stay);
nrExperiments = numel(meta);
if p.Results.newOnly
    stay = true(1,nrExperiments);
    for i=1:numel(meta)
        stay(i) = ~exists(ns.Experiment & meta(i));
    end
    meta = meta(stay);
    fullName=fullName(stay);
    nrExperiments = numel(meta);    
end
if nrExperiments ==0
    fprintf('No Neurostim files with matching paradigm and subject found in %s\n',srcFolder);
    return;
else
    fprintf('Found %d matching Neurostim files \n',nrExperiments)
    
end

if  p.Results.readFileContents
     % Read meta data from the content of the Neurostim files 
     % tmp has the meta data, c the CIC objects which are passed to
     % nsAddToDataJoint below
    for i=1:nrExperiments
        [tmp,c(i)] = ns.Experiment.readCicContents(meta(i),'root',p.Results.root);    %#ok<AGROW>         
        tmpMeta(i) = mergestruct(meta(i),tmp); %#ok<AGROW> % Merge to keep json/provenance meta.   
    end
    meta  =tmpMeta; 
else    
    c= [];
end

if p.Results.readJson
    %% Read Meta information definitions and data from JSON files.
    % Read the current meta definition for experiments
    definitionFile = fullfile(p.Results.root,['experiment_definition' metaDefinitionTag '.json']);
    if exist(definitionFile,'file')
        json = readJson(definitionFile);
        isLocked.experiment = json.Locked =="1";
        experimentMetaFields = fieldnames(json.Fields);
        if any(ismember({'starttime'},experimentMetaFields))
            error('The experiment meta definition file %s should only define new meta properties, not ''starttime''',definitionFile);
        end
        for i=1:numel(experimentMetaFields)
            [meta.(experimentMetaFields{i})] = deal("");
        end
    elseif isempty(metaDefinitionTag)
        % Use no meta definition
        isLocked.experiment = false;
    else
        error('Meta definition file %s not found',definitionFile)
    end

    % Read associated JSON files (if they exist)
    for i=1:nrExperiments
        jsonFile= regexprep(fullName{i},'(\.mat$)','.json'); % Swap extension
        if exist(jsonFile,'file')
            thisJson = readJson(jsonFile);
        else
            thisJson = struct;
        end
        metaFieldsFromJson = fieldnames(thisJson);
        for j=1:numel(metaFieldsFromJson)
            meta(i).(metaFieldsFromJson{j}) = thisJson.(metaFieldsFromJson{j});
        end
    end
end


% Convert char and cellstr to string
tExperiment = struct2table(meta,'AsArray',true);
for i=1:numel(tExperiment.Properties.VariableNames)
    if iscellstr(tExperiment.(tExperiment.Properties.VariableNames{i})) || ischar(tExperiment.(tExperiment.Properties.VariableNames{i})) %#ok<ISCLSTR>
        tExperiment.(tExperiment.Properties.VariableNames{i}) = string(tExperiment.(tExperiment.Properties.VariableNames{i}));
    end
end
% Store folder root and date in the properties so that we
% can easily reconstruct filenames later
tExperiment =addprop(tExperiment,{'root'},{'table'}) ;
tExperiment.Properties.CustomProperties.root = p.Results.root;
if p.Results.readFileContents
    tExperiment = addvars(tExperiment,c','NewVariableNames','cic'); 
end
% Sort to get a consistent order
tExperiment =  sortrows(tExperiment,{'starttime','subject','paradigm'},'ascend');
% Assign a unique ID to each row, so that we can match up with the original values after the user resorts.
tExperiment = addvars(tExperiment,strcat('#',string((1:nrExperiments)')),'NewVariableNames','id','before',1);
tExperiment = movevars(tExperiment,{'session_date','file','starttime','bytes','subject','paradigm'},'After',1);

%% B. Session Meta Data
% Retrieve definition and initialize table
tSession = unique(tExperiment(:,{'session_date','subject'}));
nrSessions = height(tSession);  
if p.Results.readJson
    definitionFile = fullfile(p.Results.root,['session_definition' metaDefinitionTag '.json']);
    if exist(definitionFile,'file')
        json = readJson(definitionFile);
        isLocked.Session = json.Locked =="1"; % The _definition file states that there should be no other meta fields.
        sessionMetaFields = fieldnames(json.Fields);
         if any(ismember({'subject','session_date'},sessionMetaFields))
            error('The session meta definition file %s should only define new meta properties, not ''session_date''',definitionFile);
        end
        nrSessionMetaFields= numel(sessionMetaFields);
    elseif isempty(metaDefinitionTag)
        nrSessionMetaFields= 0;
    else
        error('Meta definition file %s not found',definitionFile)    
    end
    emptyInit = repmat("",[height(tSession) 1]);
    if nrSessionMetaFields>0
        emtpyInitFromDef = repmat({emptyInit},[1 nrSessionMetaFields]);
        tSession = addvars(tSession,emtpyInitFromDef{:},'NewVariableNames',sessionMetaFields);
    end

    % See if there are any Session JSON files and put information in the table
    for i=1:nrSessions        
        if nrSessions==1
            sessionJsonFile = fullfile(p.Results.root,strrep(tSession.session_date,'-',filesep),tSession.subject + ".json");
             % Check to see if there are any .txt session notes (subjectNr.txt)
            txtFile = fullfile(p.Results.root,strrep(tSession.session_date,'-',filesep),tSession.subject + ".txt");
        else
            sessionJsonFile = fullfile(p.Results.root,strrep(tSession.session_date{i},'-',filesep),[tSession.subject{i} '.json']);
            % Check to see if there are any .txt session notes (subjectNr.txt)
            txtFile = fullfile(p.Results.root,strrep(tSession.session_date{i},'-',filesep),[tSession.subject{i}  '.txt']);       
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
nrSubjects = height(tSubject);
if p.Results.readJson
    definitionFile = fullfile(p.Results.root,['subject_definition' metaDefinitionTag '.json']);
    if exist(definitionFile,'file')
        json = readJson(definitionFile);
        isLocked.subject = json.Locked =="1"; % No new meta fields allowed according to _definition
        subjectMetaFields = fieldnames(json.Fields);
        if ismember('subject',subjectMetaFields)
            error('The subject meta definition file %s should only define new meta properties, not ''subject''',definitionFile);
        end
        nrSubjectMetaFields= numel(subjectMetaFields);
    elseif isempty(metaDefinitionTag)
        nrSubjectMetaFields= 0;
    else
        error('Meta definition file %s not found',definitionFile)    
    end
    emptyInit = repmat("",[height(tSubject) 1]);
    if nrSubjectMetaFields>0
        emtpyInitFromDef = repmat({emptyInit},[1 nrSubjectMetaFields]);
        tSubject = addvars(tSubject,emtpyInitFromDef{:},'NewVariableNames',subjectMetaFields);
    end
    metaDataFile =  fullfile(p.Results.root,'subject.json');
    if exist(metaDataFile,'file')
        thisJson = readJson(metaDataFile);
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
% Sort and assign a unique ID to each row.
tSubject= sortrows(tSubject,{'subject'},'ascend');
tSubject = addvars(tSubject,strcat('#',string((1:nrSubjects)')),'NewVariableNames','id','before',1);
tSubject = movevars(tSubject,{'subject'},'after',1);

if p.Results.addToDataJoint
    if p.Results.readFileContents
        cic = tExperiment.cic;
        tExperiment = removevars(tExperiment,'cic');
    else
        cic =[];
    end
    nsAddToDataJoint(tSubject,tSession,tExperiment,'cic',cic,'safeMode',p.Results.safeMode, ...
                            'root',p.Results.root,'populateFile',p.Results.populateFile,'populateMovie',p.Results.populateMovie);
end
