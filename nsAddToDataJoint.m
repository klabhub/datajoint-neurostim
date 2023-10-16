function nsAddToDataJoint(tSubject,tSession ,tExperiment,varargin)
% Takes the output of nsScan and adds new information to the DJ server.
% Depending on how nsScan was called this can add meta-data about the experiments
% only, include meta data from cic, or include all data from all plugins.
% This function is not meant to be called directly; use nsScan instead.
% INPUT
% tSubject - Output table of nsScan
% tSession - Output table of nsScan
% tExperiment - Output table of nsScan
% Parm/Value pairs
% cic - A vector of cic objects corresponding to the rows in the
% tExperiment table.
% safeMode  -   Set to false to remove tuples from the Datajoint database
%               without asking for confirmation (i.e., when updating information).
% root      - Root folder of files, defaults to getenv('NS_ROOT')
% populateFile - Run populate(ns.File) to add dependent files.
% dryrun      - Set to true to simulate what would happen without making
%               changes to the database.
%
% BK -Jan 2023

p=inputParser;
p.addParameter('safeMode',true,@islogical);
p.addParameter('root',getenv('NS_ROOT'));
p.addParameter('populateFile',true,@islogical)
p.addParameter('populateMovie',true,@islogical)
p.addParameter('dryrun',false,@islogical)
p.addParameter('cic',[]); % A vector of cic objects. One per row of the experiment table.
p.parse(varargin{:});

currentSafeMode= dj.config('safemode');
dj.config('safemode',p.Results.safeMode);
warnstate =warning('query');
warning('off', 'DataJoint:longCondition');

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
[newExpts,~,newTplsRows] = insertNewTuples(tExperiment,ns.Experiment,p.Results.dryrun);
if  ~p.Results.dryrun && ~isempty(newExpts) 
    if ~isempty(p.Results.cic)
        updateWithFileContents(ns.Experiment & newExpts,p.Results.cic(newTplsRows))
    end
    if p.Results.populateFile
        % Populate the File table (all files associated with this experiment)
        populate(ns.File, newExpts);
    end
    if p.Results.populateMovie
        % Populate the Movie table (subset of files with a specific set of extensions known to be movies) 
        populate(ns.Movie,newExpts);
    end
end

% Restore setting
dj.config('safemode',currentSafeMode);
warning(warnstate);


end
function [newTpls,newMetaTpls,newTplsRows] = insertNewTuples(tbl,djTbl,dryrun)
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
newMetaTpls = [];  %
updateMetaTpls = []; % With an empty metaTpl defined we don't have to check for nr>0 below

pkey = djTbl.primaryKey;
hdr  = djTbl.header;
djTblName = djTbl.className;
state= warning('query');
tblFields= intersect(tbl.Properties.VariableNames,hdr.names); % Columns in the DJ table
metaFields = setdiff(tbl.Properties.VariableNames,hdr.names); % Columns not in the DJ table
% Convert all meta and dj variables that are defined string to string in
% the tbl.
warning('off','MATLAB:table:convertvars:ConvertCharWarning')
tblFieldAttributes = [hdr.attributes];
tblFieldAttributes = tblFieldAttributes(ismember({tblFieldAttributes.name},tblFields));
tblFieldName = {tblFieldAttributes.name};
stringFields = cat(2,metaFields, tblFieldName([tblFieldAttributes.isString]));
tbl = convertvars(tbl,stringFields,'string');
warning(state);

nrMeta = numel(metaFields);
% Hack, find the associated meta.SubejctMeta, ExperimentMeta,etc.
djMetaTblName  = [djTbl.className 'Meta'];
djMetaTbl = feval(djMetaTblName);


% Loop over the new table
newTplsRows = [];
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
        newTplsRows = cat(1,newTplsRows,row);
        if nrMeta>0
            newMetaTpls= cat(1,newMetaTpls,thisMetaTpl);       
        end
    else
        % Existing tuple.
        % Check whether the new one is different from the
        % existing one
        if exists(djTbl & table2struct(tbl(row,tblFields))) && nrMeta>0
            % Tpls are the same, check if the meta information is different.
            fromDbase = ns.getMeta(djTbl & dbTpl,metaFields);            
            fromDbase = convertvars(fromDbase,1:width(fromDbase),"string"); % Match the "" format of the tbl
            if isempty(fromDbase) 
                % Everything is new
                newMetaTpls= cat(1,newMetaTpls,thisMetaTpl);
            elseif  ~isempty(setdiff(cat(2,pkey,metaFields),fromDbase.Properties.VariableNames)) 
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
            newTplsRows = cat(1,newTplsRows,row);
            % The update will remove meta data, so add those as new as well
            if nrMeta>0
                newMetaTpls= cat(1,newMetaTpls,thisMetaTpl);
            end
        end
    end
end
% Remove empty meta data to avoid filling the database with empty entries.
% This removal is only done for new meta tuples; updates can be empty to
% remove the old information.
if ~isempty(newMetaTpls)
    stay = ~cellfun(@isempty,{newMetaTpls.meta_value},'uni',true);
    newMetaTpls = newMetaTpls(stay);
end


if dryrun
    fprintf('[DRYRUN] %s: %d new , %d updated \n %s: %d new meta, %d meta updated\n',djTblName,numel(newTpls),numel(updateTpls),djMetaTblName,numel(newMetaTpls),numel(updateMetaTpls))
else
    % Would be nice to wrap in a transactiomn, but cannot insert before
    % commiting the del.
    try
        fprintf('Updating DataJoint for %s ...\n',djTblName)
        if ~isempty(newTpls)
            fprintf('Adding new tuples to %s \n',djTblName)
            insert(djTbl,newTpls)
            fprintf('\t Done. %d new tuples \n',numel(newTpls))
        end

        if ~isempty(updateTpls)
            fprintf('Updating tuples in %s \n',djTblName)
            if true
                % Update one tpl and one field at a time.
                % In theory this could affect referential integrity, but the
                % delete/replace approach is too cumbersome (e.g., renaming one
                % subject would delete all data associated with that subject)
                fieldsToUpdate = setdiff(fieldnames(updateTpls),pkey)';
                for tpl =updateTpls
                    thisDj= djTbl & ns.stripToPrimary(djTbl,tpl);
                    for fld = fieldsToUpdate
                        update(thisDj,fld{1},tpl.(fld{1}))
                    end
                end
            else
                del(djTbl & ns.stripToPrimary(djTbl,updateTpls)); % Deletes tpl and associated meta
                insert(djTbl,updateTpls)
            end
            fprintf('\t Done. %d updated tuples\n',numel(updateTpls))
        end

        if ~isempty(newMetaTpls)
            fprintf('Adding new Meta tuples in %s \n',djMetaTblName)
            insert(djMetaTbl,newMetaTpls)
            fprintf('\t Done. %d new tuples\n',numel(newMetaTpls))
        end
        if ~isempty(updateMetaTpls)
            fprintf('Updating Meta tuples in %s \n',djMetaTblName)
            % Because the collection of meta data could have changed, we delete
            % all, then replace to get the current set of meta fields and
            % values.
            del(djMetaTbl & ns.stripToPrimary(djMetaTbl,updateMetaTpls)); % Delete meta
            insert(djMetaTbl,updateMetaTpls); % Insert
            fprintf('\t Done. %d updated tuples\n',numel(updateMetaTpls))
        end
    catch me
        fprintf(2,'Failed to insert:')
        tbl
        me.message
    end
end
newTpls = cat(1,newTpls,updateTpls); % Both are considered new as they may need postprocessing(e.g. for Experiment table; read the file or populate ns.File).
newMetaTpls= cat(1,newMetaTpls,updateMetaTpls);
end


