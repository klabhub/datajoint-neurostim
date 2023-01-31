function nsAddToDataJoint(tSubject,tSession ,tExperiment,varargin)
% Takes the output of nsScan and adds new information to the DJ server.
% INPUT
% tSubject - Output table of nsScan
% tSession - Output table of nsScan
% tExperiment - Output table of nsScan
% Parm/Value pairs
% 'readFileContents' - Set to true to read each experiment file when adding
%                       it to the database (To add information obtained
%                       from its content). [false]
% safeMode  -   Set to false to remove tuples from the Datajoint database
%               without asking for confirmation (i.e., when updating information).
% root      - Root folder of files, 
% BK -Jan 2023

p=inputParser;
p.addParameter('readFileContents',false,@islogical);
p.addParameter('safeMode',true,@islogical);
p.addParameter('root',get(ns.Global,'root'));
p.parse(varargin{:});

currentSafeMode= dj.config('safemode');
dj.config('safemode',p.Results.safeMode);


%% Add Subjects 
% Replace emppty dates with 0 date
if ismember('dob',tSubject.Properties.VariableNames)
    tSubject{tSubject.dob=="",'dob'} = "0000-00-00";
end
tSubject= removevars(tSubject,'id');
insertNewTuples(tSubject,ns.Subject);

%% Add Sessions

tSession= removevars(tSession,'id');
insertNewTuples(tSession,ns.Session);

%% Add Experiments
tExperiment = removevars(tExperiment,'id');
[newExpts] = insertNewTuples(tExperiment,ns.Experiment);
if p.Results.readFileContents
    % Will read each (new) file and add its contents to DataJoint
    updateWithFileContents(ns.Experiment & newExpts,'root',tExperiment.Properties.CustomProperties.root);
end

% Restore setting
dj.config('safemode',currentSafeMode);



end
function [tpl,metaTpl] = insertNewTuples(tbl,djTbl)
% Given a table read from a folder and a table (Relvar) from the Datajoint
% databse, determin which tuples are new and the meta data associated with
% those new tuples and insert those in the Datajoint tables.
%
% INPUT 
% tbl - Table with information that nsScan creates from filenames plus json files
%               This is a Matlab table. 
% djTbl - dj.Relvar with the information already in the database.
% 
% OUTPUT
% tpl  - Tuples with new information
% metaTpl - Tuples with new meta information.

tpl= [];
metaTpl= [];
pkey = djTbl.primaryKey;
inFiles = unique(tbl(:,pkey));
inDatabase = fetchtable(djTbl);
% Determine which rows to keep
if ~isempty(inDatabase)
    [~,ix] = setdiff(inFiles,inDatabase);
else % All 
    ix = 1:height(inFiles);
end
if isempty(ix)
    fprintf('No new tuples for %s \n',djTbl.className)
    
    return;
end
% Determine which columns to keep
hdr = djTbl.header;
[overlapFields,ixKeep]  = intersect(tbl.Properties.VariableNames,hdr.names);
%convertvars(tbl(ix,overlapFields),tbl.Properties.VariableNames(ixKeep),'char')
tpl = table2struct(tbl(ix,overlapFields)); % Only the fields that are defined in the SQL database

% The cols that are not defined in SQL can be added to a linked meta data table
isMeta = setdiff(tbl.Properties.VariableNames,hdr.names);
metaTbl = table;
for i=ix(:)'
    for j= i:numel(isMeta)
        metaTbl = [metaTbl; [tbl(i,pkey) table(isMeta(j), tbl{i,isMeta{j}},'VariableNames',{'meta_name','meta_value'})]]; %#ok<AGROW> 
    end
end
%convertvars(metaTbl,metaTbl.Properties.VariableNames,'char')
metaTpl = table2struct(metaTbl);
fprintf('Adding new tuples to %s \n',djTbl.className)
insert(djTbl,tpl)
fprintf('\t Done. %d new tuples\n',numel(tpl),djTbl.className)
if ~isempty(metaTpl)
    % Hack, find the associated meta.SubejctMeta, ExperimentMeta,etc.
    metaTblName  = [djTbl.className 'Meta'];
    metaTable = feval(metaTblName);
    insert(metaTable,metaTpl);
    fprintf('Adding %d new tuples to %s \n',numel(metaTpl),metaTblName);
end
end