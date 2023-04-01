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
% root      - Root folder of files, defaults to getenv('NS_ROOT')
% populateFile - Run populate(ns.File) to add dependent files.
% cicOnly   - Set tot true to only add CIC information (and not other
%               plugins)
% dryrun      - Set to true to simulate what would happen without making
%               changes to the database.
% BK -Jan 2023

p=inputParser;
p.addParameter('readFileContents',false,@islogical);
p.addParameter('safeMode',true,@islogical);
p.addParameter('root',getenv('NS_ROOT'));
p.addParameter('populateFile',true,@islogical)
p.addParameter('dryrun',false,@islogical)
p.addParameter('cicOnly',false,@islogical)
p.parse(varargin{:});

currentSafeMode= dj.config('safemode');
dj.config('safemode',p.Results.safeMode);


%% Add Subjects
% Replace emppty dates with 0 date
if ismember('dob',tSubject.Properties.VariableNames)
    tSubject{tSubject.dob=="",'dob'} = "0000-00-00";
end
tSubject= removevars(tSubject,'id');
insertNewTuples(tSubject,ns.Subject,p.Results.dryrun);

%% Add Sessions

tSession= removevars(tSession,'id');
insertNewTuples(tSession,ns.Session,p.Results.dryrun);

%% Add Experiments
tExperiment = removevars(tExperiment,'id');
[newExpts] = insertNewTuples(tExperiment,ns.Experiment,p.Results.dryrun);
if  ~p.Results.dryrun && ~isempty(newExpts)
    if p.Results.readFileContents
        % Read each (new) file and add its contents to DataJoint
        % In this process the Experiment tpl is first deleted then re-added
        % with the information from the file. (No code yet to avoid this
        % add/delete/add)
        updateWithFileContents(ns.Experiment & newExpts,'root',tExperiment.Properties.CustomProperties.root,'cicOnly',p.Results.cicOnly);
    end

    if p.Results.populateFile
        populate(ns.File, newExpts)
    end
end

% Restore setting
dj.config('safemode',currentSafeMode);


end
function [newTpls,newMetaTpls] = insertNewTuples(tbl,djTbl,dryrun)
% Given a table read from a folder and a table (Relvar) from the Datajoint
% databse, determin which tuples are new and the meta data associated with
% those new tuples and insert those in the Datajoint tables.
%
% INPUT
% tbl - Table with information that nsScan creates from filenames plus json files
%               This is a Matlab table.
% djTbl - dj.Relvar with the information already in the database.
% dryrun - Set to true to do a dryrun
% OUTPUT
% tpl  - Tuples with new information
% metaTpl - Tuples with new meta information.

%% Determine whether there are entries that are not yet in the database
newTpls = [];
updateTpls = [];
newMetaTpls = [];
updateMetaTpls = [];

pkey = djTbl.primaryKey;
hdr  = djTbl.header;
djTblName = djTbl.className;
state= warning('query');
warning('off','MATLAB:table:convertvars:ConvertCharWarning')
tbl = convertvars(tbl,1:width(tbl),'string');
warning(state);

tblFields= intersect(tbl.Properties.VariableNames,hdr.names); % Columns in the DJ table
metaFields = setdiff(tbl.Properties.VariableNames,hdr.names); % Columns not in the DJ table
nrMeta = numel(metaFields);
% Hack, find the associated meta.SubejctMeta, ExperimentMeta,etc.
djMetaTblName  = [djTbl.className 'Meta'];
djMetaTbl = feval(djMetaTblName);


% Loop over the new table
for row=1:height(tbl)
    thisPrimaryTpl =table2struct(tbl(row,pkey));
    dbTpl = fetch(djTbl & thisPrimaryTpl);
    thisTblTpl =table2struct(tbl(row,tblFields));% Only the fields that are defined in the SQL database.
    for m=1:nrMeta
        thisMetaTpl(m,1) = mergestruct(thisPrimaryTpl,struct('meta_name',metaFields{m},'meta_value',char(tbl{row,metaFields{m}}))); %#ok<AGROW>
    end
    if isempty(dbTpl)
        % No match with the primary key
        % Add to the newTpls array
        newTpls = cat(1,newTpls,thisTblTpl);
        newMetaTpls= cat(1,newMetaTpls,thisMetaTpl);
    else
        % Existing tuple.
        % Check whether the new one is different from the
        % existing one
        if exists(djTbl & table2struct(tbl(row,tblFields)))
            % Tpls are the same, check if the meta information is different.
            fromDbase = ns.getMeta(djTbl & dbTpl,metaFields);
            fromDbase = convertvars(fromDbase,1:width(fromDbase),"string"); % Match the "" format of the tbl
            if isempty(fromDbase)
                % Everything is new
                newMetaTpls= cat(1,newMetaTpls,thisMetaTpl);
            elseif  ~isequal(fromDbase.Properties.VariableNames,cat(2,pkey,metaFields))
                % Some newly defined meta fields. Update all
                updateMetaTpls= cat(1,updateMetaTpls,thisMetaTpl);
            else 
                % Same meta fields, potentially different values
                [~,ix] = setdiff(tbl(row,cat(2,pkey,metaFields)),fromDbase);
                updateMetaTpls= cat(1,updateMetaTpls,thisMetaTpl(ix));
            end
        else
            %Add to the updateTpls array
            updateTpls = cat(1,updateTpls, thisTblTpl);
            % The update will remove meta data, so add those as new as well
            newMetaTpls= cat(1,newMetaTpls,thisMetaTpl);
        end
    end
end

if dryrun
    fprintf('[DRYRUN] %s: %d new , %d updated \n %s: %d new meta, %d meta updated\n',djTblName,numel(newTpls),numel(updateTpls),djMetaTblName,numel(newMetaTpls),numel(updateMetaTpls))
else
    % Would be nice to wrap in a transactiomn, but cannot insert before
    % commiting the del.

    fprintf('Updating DataJoint for %s ...\n',djTblName)
    if ~isempty(newTpls)
        fprintf('Adding new tuples to %s \n',djTblName)
        insert(djTbl,newTpls)
        fprintf('\t Done. %d new tuples \n',numel(newTpls))
    end

    if ~isempty(updateTpls)
        fprintf('Updating tuples in %s \n',djTblName)
        del(djTbl & ns.stripToPrimary(djTbl,updateTpls)); % Deletes tpl and associated meta
        insert(djTbl,updateTpls)
        fprintf('\t Done. %d updated tuples\n',numel(updateTpls))
    end

    if ~isempty(newMetaTpls)
        fprintf('Adding new Meta tuples in %s \n',djMetaTblName)
        insert(djMetaTbl,newMetaTpls)
        fprintf('\t Done. %d new tuples\n',numel(newMetaTpls))
    end
    if ~isempty(updateMetaTpls)
        fprintf('Updating Meta tuples in %s \n',djMetaTblName)
        del(djMetaTbl & ns.stripToPrimary(djMetaTbl,updateMetaTpls)); % Delete meta
        insert(djMetaTbl,updateMetaTpls); % Insert
        fprintf('\t Done. %d updated tuples\n',numel(updateMetaTpls))
    end
end
newTpls = cat(1,newTpls,updateTpls); % Both are considered new as they may need postprocessing(e.g. for Experiment table; read the file or populate ns.File).
newMetaTpls= cat(1,newMetaTpls,updateMetaTpls);
end


