function setMeta(tbl,metaParm,metaVal,comments,pv)
% Function to set or add meta data to json files corresponding to
% experiments, sessions, or subjects.
arguments
    tbl
    metaParm (1,:) string
    metaVal (1,:) string
    comments (1,1) string  = ""
    pv.dryrun (1,1) logical = true
end

% Sanity checks and singleton expansion
nrTpls = count(tbl);
assert(isscalar(metaParm) || nrTpls==numel(metaParm),'Number of rows in the table must match the number of meta parameters');
assert(isscalar(metaVal) || nrTpls==numel(metaVal), 'Number of rows in the table must match the number of meta parameters');
if isscalar(metaParm)
    metaParm = repmat(metaParm,[nrTpls 1]);
end
if isscalar(metaVal)
    metaVal= repmat(metaVal,[nrTpls 1]);
end

% Loop over the table
tplCntr = 0;
if pv.dryrun
    fprintf('DRYRUN:\n');
end

for tpl = fetch(tbl)'
    tplCntr = tplCntr +1;
    % Determine json file
    switch class(tbl)
        case "ns.Subject"
            error('niy')
        case "ns.Experiment"            
            [fldr,thisFile,~] = fileparts(file(ns.Experiment &tpl));            
        case "ns.Session"
            error('niy')
        otherwise
            error('Meta data for %s???',class(tbl));
    end
    jsonFile =fullfile(fldr ,thisFile +".json");
    % Read or create metaData
    if exist(jsonFile,"file")
        metaData = readJson(jsonFile);
        changed  = ~isfield(metaData,metaParm(tplCntr)) || metaVal(tplCntr)~= metaData.(metaParm(tplCntr));
        if changed
            metaData(1).(metaParm(tplCntr)) = metaVal(tplCntr);
            if comments ~=""
                if isfield(metaData(1),"comments") && metaData(1).comments ~=""                     
                    % Append with ;
                    metaData(1).comments =  metaData(1).comments + ";" + comments;
                else
                    % Create new or set
                    metaData(1).comments = comments;                    
                end               
            end
        end
    else
        changed = true;
        metaData = struct(metaParm(tplCntr),metaVal(tplCntr),"comments",comments);        
    end
    % Save changed meta data
    if changed
        % Write the updated JSON back to the file
        fprintf('Update meta data for %s in %s to %s\n',metaParm(tplCntr),jsonFile, metaVal(tplCntr))
        if ~pv.dryrun              
            writeJson(metaData,jsonFile);
        end
    else
        fprintf('No change in meta data for %s in %s\n',metaParm(tplCntr),jsonFile)
    end
end
if pv.dryrun
    fprintf('DRYRUN Complete.\n');
end
end