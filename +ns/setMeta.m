function setMeta(tbl,metaParm,metaVal,comments)
arguments
    tbl
    metaParm (1,:) string
    metaVal (1,:) string
    comments (1,1) string  = ""
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
for tpl = fetch(tbl)'
    tplCntr = tplCntr +1;
    % Determine json file
    switch class(tbl)
        case "ns.Subject"
            error('niy')
        case "ns.Experiment"
            jsonFile = strrep(file(ns.Experiment &tpl),'.mat','.json');
        case "ns.Session"
            error('niy')
        otherwise
            error('Meta data for %s???',class(tbl));
    end
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
        fprintf('Update meta data for %s in %s\n',metaParm(tplCntr),jsonFile)
        writeJson(metaData,jsonFile);
    else
        fprintf('No change in meta data for %s in %s\n',metaParm(tplCntr),jsonFile)
    end
end
end