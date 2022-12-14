function djScan(varargin)
% This function scans filenames in (and around) the folder correspondong
% to a certain date under a given root folder and adds Neurostim experiments
% found in those folders to the datajoint database.
%
% This works because Neurostim always (?) generates files in the following
% format:
% root\YYYY\MM\dd\subject.paradigm.startTime.mat
%
% INPUT
% 'root' - Top level folder (e.g. Z:\). Folders below this root folder should
% correspond to the years. [Read from ns.Global 'root' property by default]
% date - Date to scan.  [today]
% schedule - 'y' - Scan the year in which the date falls.
%          - 'm' - Scan the month
%          - ['d'] - Scan the specified day only.
% readFileContents - Read each Neurostim file and update the database with its
%               values [false]
%           With this set to false, this creates a quick overview of files
%           in the data root folder. File content can later be added using
%           the updateWithFileContents(ns.Experiment)
%
% ignore  - File extensions to ignore {'.ini','.cache'}
% safemode   = [true] -Ask confirmation before dropping tables in the
%                       database
% paradigms = Cell array of paradigms to include. Leave empty to include
% all. [{}]
% fileType -  File extenstions to scan [.mat]
% folderFun -  Function, provided by the user, to handle special folders.
%           The djScan function will call this function with a list of
%           folders. The user is responsible for returning experiments from
%           the relevant folders. See djScanDicomFolder for an example. Not
%           commonly needed. []. 
% OUTPUT
% None
%
% BK - April 2022.

rt = fetchn(ns.Global & 'name=''root''' ,'value','ORDER BY id DESC LIMIT 1');
if isempty(rt)
    rt = pwd;
else
    rt =rt{1};
end
p = inputParser;
p.addParameter('root',rt);
p.addParameter('date',now);
p.addParameter('schedule','d');
p.addParameter('readFileContents',false);
p.addParameter('ignore',{'.ini','.cache'});
p.addParameter('safemode',true);
p.addParameter('paradigms',{});
p.addParameter('fileType','.mat')
p.addParameter('folderFun','');
p.parse(varargin{:});

dj.config('safemode',p.Results.safemode);


switch (p.Results.schedule)
    case 'y'
        % Allow more than one folder below year
        srcFolder = fullfile(p.Results.root,datestr(p.Results.date,'yyyy'),'**');
    case 'm'
        % Exactly one folder below month
        srcFolder = fullfile(p.Results.root,datestr(p.Results.date,'yyyy/mm'),'**');
    case 'd'
        % Inside the day folder
        srcFolder = fullfile(p.Results.root,datestr(p.Results.date,'yyyy/mm/dd'));
    otherwise
        error('Unknown schedule %s',p.Results.schedule)
end

fprintf('Scanning %s ...\n',srcFolder)

% Find the files matching the wildcard
files= dir(fullfile(srcFolder,['*' p.Results.fileType]));
pathFromRoot = strrep({files.folder},p.Results.root,'');
fullName = strcat(pathFromRoot',filesep,{files.name}');
% Prune to get files that Neurostim generates
if strcmpi(filesep','\')
    fs = '\\';
else
    fs = filesep;
end
pattern = ['(?<date>\d{4,4}' fs '\d{2,2}' fs '\d{2,2})' fs '(?<subject>\w{2,10})\.(?<paradigm>\w+)\.(?<startTime>\d{6,6})\.'];
nsDataFiles = regexp(fullName,pattern,'names');
% Prune those file that did not match
out = cellfun(@isempty,nsDataFiles);
nsDataFiles= [nsDataFiles{~out}]; % Create a struct array with the relevant info
if ~isempty(nsDataFiles)
    [nsDataFiles.folder] = deal(false); % Mark as true NS files, not folders with other files.
    [nsDataFiles.type] = deal(p.Results.fileType);
end

if ~isempty(p.Results.paradigms) && ~isempty(nsDataFiles)
    out = ~ismember({nsDataFiles.paradigm},p.Results.paradigms);
    if any(out)
        fprintf('Skipping %d files with non-matching paradigms\n',sum(out))        
        nsDataFiles(out)=[];        
    end
end


if ~isempty(p.Results.folderFun)
    % Special handling of sub folders. Not usually needed, but see djScanDicomFolder
    % for an example of how to include DICOM MRI data folders.
    % Note that all files in one of these folders (i.e. those returned by
    % folderFun) will be added to the ns.Files table.
    folders= dir(fullfile(srcFolder,'*'));
    folders(~[folders.isdir])=[];  % Remove non folders  
    folderFullName = strcat({folders.folder}',filesep,{folders.name}');
    nsDataFilesFromFolder = p.Results.folderFun(folderFullName);
    if ~isempty(nsDataFilesFromFolder)
        nsDataFiles = cat(2,nsDataFiles,nsDataFilesFromFolder);
    end
end



if isempty(nsDataFiles)
    fprintf('No files found in %s\n',p.Results.root);
    return;
end

fprintf('Foound %d files with matching paradigms\n',numel(files))

%% Add the new subjects (if any)
% Assuming that all non-primary keys are nullable.
uSubjects = unique({nsDataFiles.subject});
tbl  = ns.Subject;
knownSubjects = fetch(tbl,'subject');
newSubjects = setdiff(uSubjects,{knownSubjects.subject});
nullFields = tbl.nonKeyFields;
tmp = cell(1,2*numel(nullFields));
[tmp{1:2:end}] = deal(nullFields{:});
[tmp{2:2:end}] = deal([]);
newSubjects = struct(tbl.primaryKey{1},newSubjects,tmp{:});
fprintf('Adding %d new subjects \n',numel(newSubjects))
insert(ns.Subject,newSubjects)

%% Loop over datafiles (which should correspond to Neurostim experiments or folders with other files)
nrDataFiles  = numel(nsDataFiles);
for i=1:nrDataFiles
    if nsDataFiles(i).folder
        fname = fullfile(nsDataFiles(i).date,sprintf('%s%s',nsDataFiles(i).subject,nsDataFiles(i).type),nsDataFiles(i).paradigm);
    else
        fname = fullfile(nsDataFiles(i).date,sprintf('%s.%s.%s%s',nsDataFiles(i).subject,nsDataFiles(i).paradigm,nsDataFiles(i).startTime,nsDataFiles(i).type));
    end
    fprintf('Processing file %d of %d (%s)\n',i,nrDataFiles,fname)

    %% Find or add Session
    qry =struct('session_date', datestr(nsDataFiles(i).date,29),...  % Convert to ISO 8601
        'subject',nsDataFiles(i).subject);
    thisSession = ns.Session & qry;
    if ~thisSession.exists
        insert(ns.Session,qry);
    end

    %% Find or add Experiment
    qry.starttime = datestr(datenum(nsDataFiles(i).startTime,'HHMMSS'),'HH:MM:SS');
    qry.paradigm = nsDataFiles(i).paradigm;
    thisExperiment = ns.Experiment  & qry;
    file = fullfile(p.Results.root,fname);
    if ~thisExperiment.exists || p.Results.readFileContents
        key = qry;
        key.file    = file;
        if p.Results.readFileContents && ~exist(file,'dir')
            updateWithFileContents(ns.Experiment,key);
        else
            % Just insert the file info
            try
                insert(ns.Experiment,key);
            catch
                fprintf(2,"Duplicate entry? Starttime incorrect ?\n")
            end
        end
    end
    qry = rmfield(qry,'paradigm');

    %% Find or add Files  - this requires a search for all files that have the same prefix

    if nsDataFiles(i).folder
        % This is an entry that represents a folder from
        % p.Results.folderFun, not a neurostim experiment .mat file.  Add
        % all files in the folder
         linkedFiles= dir(fullfile(file,'*'));
         pth = file;
    else
        % This is a neurostim experiment. Search for files with matching prefix (subject.paradigm.startTime.*) in 
        % the same folder and all files in the folder with this name
        [pth,f] =fileparts(file);
        prefix = fullfile(pth,[f '*']);
        inFolder = dir(prefix);
        prefix = fullfile(prefix,'*');
        inSubFolder = dir(prefix);
        linkedFiles = cat(1,inFolder,inSubFolder);
    end
    
    % Remove folders
    linkedFiles([linkedFiles.isdir]) =[];
    % Remove common junk files
    match = regexp({linkedFiles.name},'.*\.(?<ext>.*$)','names');
    match = [match{:}];
    ext = strcat('.',{match.ext});
    out = ismember(ext,p.Results.ignore);
    linkedFiles(out) =[];
    ext(out) = [];
    
    % Add each one
    for f=1:numel(linkedFiles)
        relFolder = strrep(linkedFiles(f).folder,pth,'');     
        qry.filename = fullfile(relFolder,linkedFiles(f).name);
        thisFile = ns.File & qry;
        if ~thisFile.exists
            qry.extension = ext{f};
            insert(ns.File,qry);
            qry = rmfield(qry,"extension");
        end
    end
end

% Show an overview of the current database
djReport;
end








