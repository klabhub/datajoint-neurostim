function scanMeta(tbl,pv)
% Scan the file tree for JSON files with meta data and update the datajoint
% database with the new values.
%
% Typically these data will have been entered in the nsMeta app (by hand),
% this funciton synchronizes the state of the JSON files in the file tree
% with the tables in datajoint. Because this only runs over
% subjects/session/experiments that are already in the database, it will
% not add new tpls.
%
% INPUT
% tbl  -  A subject,session, or experiment table
arguments
    tbl (1,1) {mustBeA(tbl,{'ns.Experiment','ns.Subject','ns.Session'})}
    pv.dryrun (1,1) logical = false
    pv.NS_ROOT = getenv("NS_ROOT")
end
% Get the meta table
className = class(tbl);
metaTbl = feval([className 'Meta']) & tbl;

%% Special handling of ns.Subject
% dob,sex,and species are stored in the main table (not SubjectMeta).
% Handle them as updates here (slow).
if strcmpi(className ,'ns.Subject')
    % Read the single json
    jsonFile  =fullfile(pv.NS_ROOT,'subject.json');
    allJson = readJson(jsonFile);
    % Limit to subjects already in the database with different values
    % But this only checks for "meta" datain the ns.Subject table ; the
    % code below runs over all the json meta data without this restriction
    needsUpdate = ns.Subject - allJson;
    needsUpdate = fetchtable(needsUpdate);
    needsUpdate = ismember([allJson.subject],needsUpdate.subject);
    jsonNeedsUpdate  =allJson(needsUpdate);    
    fprintf('Updating %d subject entries.',numel(jsonNeedsUpdate))
    for j =1:numel(jsonNeedsUpdate)
        key = tbl & struct('subject',jsonNeedsUpdate(j).subject);      
        if jsonNeedsUpdate(j).dob==""
            update(key,'dob',[]);
        else
            update(key,'dob',char(jsonNeedsUpdate(j).dob));
        end
        if jsonNeedsUpdate(j).sex==""
            update(key,'sex','u');
        else
            update(key,'sex',jsonNeedsUpdate(j).sex);
        end
        update(key,'species',char(jsonNeedsUpdate(j).species));
        fprintf('.');
        if mod(j,40)==0;fprintf('\n');end
    end
    fprintf('Done.\n');
end

%% Handle meta data stored in the meta table
% (For Session and Experiment these are the only ones)
fprintf('Constructing meta data to replace. \n')
dj.conn().startTransaction
try
    if pv.dryrun
        fprintf('DRY RUN - NO CHANGES WILL BE MADE\n')
    else
        % Delete the old meta table - everything is stored in JSON which will
        % be read below
        delQuick(metaTbl &tbl);
    end
    % Loop over the table to read the json files and construct an array of tpls
    for key  = fetch(tbl)'
        switch (className)
            case 'ns.Subject'
                thisJson = allJson([allJson.subject]==key.subject);
                if isempty(thisJson)
                    fprintf("No metadata found in %s for %s \n ",jsonFile, key.subject)
                else
                    thisJson = rmfield(thisJson,{'subject','dob','species','sex'}); % These fields are stored in the ns.Subject table and have been updated above
                end
            case 'ns.Session'
                jsonFile  =fullfile(folder(ns.Session&key),[key.subject '.json']);
                thisJson = readIfExists(jsonFile);
            case 'ns.Experiment'
                [fldr,matFile]  =fileparts(file(ns.Experiment & key));
                jsonFile= fullfile(fldr,matFile + ".json");
                thisJson = readIfExists(jsonFile);
        end
        if ~isempty(thisJson)
            fn=fieldnames(thisJson);
            newTpl = mergestruct(key,struct('meta_name',fn,...
                'meta_value',struct2cell(thisJson)));
            if exist('tpl',"var")
                tpl = catstruct(1,tpl,newTpl);
            else
                tpl = newTpl;
            end
        end
    end
    % Add the found tpls to the data base.
    if exist('tpl',"var")
        fprintf('Inserting %d meta data.\n' ,numel(tpl))
        if pv.dryrun
            fprintf('DRY RUN - NO CHANGES MADE.\n')
        else
            % In with the new
            insert(metaTbl,tpl)
        end
    else
        fprintf('No relevant meta data found. \n')
    end

catch me
    fprintf(2,'An error occurred. Cancelling the transaction. No changes have been made. (%s)\n',me.message)
    dj.conn().cancelTransaction
    return;
end
dj.conn().commitTransaction
end

%% Helper
function thisJson = readIfExists(f)
if exist(f,"file")
    thisJson = readJson(f);
else
    fprintf("No JSON metadata found. %s does not exist\n ",f)
    thisJson = [];
end
end




