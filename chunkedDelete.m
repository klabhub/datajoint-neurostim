function chunkedDelete(targetQuery, batchSize, extraRelvars, autoDetectFK)
% To avoid long server locks, delete rows in tables in batches.
% This will first scan the query tree, and then delete rows from the bottom
% up to avoid long table locks.
% If dj.congif('safeMode') =1 the function asks for confirmation (once, after showing
% the deletion plan)/
arguments
    targetQuery (1,1) dj.internal.GeneralRelvar
    batchSize (1,1) double = 500 %Delete this many rows each call.
    extraRelvars cell = {} % Optional cell array of additional relvars (possibly cross-schema) to include.
    autoDetectFK (1,1) logical = true % When true, include external FK children discovered via INFORMATION_SCHEMA
end

if ~exists(targetQuery)
    fprintf('Nothing to delete. \n');
    return;
end

% 1.Analysis
disp('Analyzing dependency tree and tuple counts...');

% Collect all descendants, including optional cross-schema relvars.
descBase = targetQuery.descendants;
descExtra = {};
for k = 1:numel(extraRelvars)
    rv = extraRelvars{k};
    if ischar(rv) || isstring(rv)
        rv = feval(rv);
    end
    descExtra = [descExtra, rv.descendants]; %#ok<AGROW>
end

list = unique([descBase, descExtra], 'stable');

% Optionally pull in external FK children by inspecting INFORMATION_SCHEMA.
if autoDetectFK
    listFullNames = cellfun(@(c) feval(c).fullTableName, list, 'uni', false);
    extChildren = discoverExternalChildren(listFullNames);
    list = unique([list, extChildren], 'stable');
end

list = flip(list);


totalTuples = 0;
summary = cell(length(list), 2);
for i = 1:length(list)
    tableObj = feval(list{i});
    toDelete = tableObj & targetQuery.proj();
    c = count(toDelete);
    summary{i, 1} = list{i};
    summary{i, 2} = c;
    totalTuples = totalTuples + c;
end


% Display Plan
fprintf('\n--- DELETE PLAN (Bottom-Up) ---\n');
for i = 1:size(summary, 1)
    if summary{i, 2} > 0
        fprintf('%8d tuples from %s\n', summary{i, 2}, summary{i, 1});
    end
end

% 2. Safety Mode Confirmation
if dj.config('safemode') && totalTuples > 0
    prompt = sprintf('\nSafemode is ON. Delete %d tuples? (y/n): ', totalTuples);
    if ~strcmpi(input(prompt, 's'), 'y')
        disp('Aborted.'); return;
    end
end

% 3. Resilient Execution
for i = 1:length(list)
    tableName = list{i};
    tableObj = feval(tableName);
    toDelete = tableObj & targetQuery.proj();
    numTuples = count(toDelete);
    if numTuples == 0; continue; end

    fprintf('\nProcessing %s...\n', tableName);
    keys = fetch(toDelete); % Fetch keys into MATLAB RAM

    j = 1;
    while j <= numTuples
        endIdx = min(j + batchSize - 1, numTuples);
        currentBatch = keys(j:endIdx);

        try
            delQuick(tableObj & currentBatch); % Commit-per-batch
            fprintf('  Deleted %d/%d\n', endIdx, numTuples);
            j = j + batchSize; % Advance only on success
        catch ME
            if contains(ME.message, 'gone away') || contains(ME.message, 'Lost connection')
                fprintf('\n[CONNECTION LOST] Attempting to reconnect...\n');
                pause(5); % Give the server a moment to finish restarting
                try
                    dj.conn.reopen(); % Re-establish the pipe
                    fprintf('[RECONNECTED] Resuming batch at index %d...\n', j);
                    % Loop repeats with the same j index
                catch
                    error('Fatal: Could not reconnect to MariaDB.');
                end
            else
                rethrow(ME); % Crash on actual SQL syntax errors
            end
        end
    end
end
disp('Resilient bottom-up deletion complete.');
end


function extList = discoverExternalChildren(fullTableNames)
% Query INFORMATION_SCHEMA to find tables (any schema) that reference the
% provided tables, then return their class names (including their own
% descendants) so they can be deleted first.

    extList = {};
    if isempty(fullTableNames); return; end

    % Build an IN (...) list like ('`db`.`table`','`db2`.`table2`')
    quoted = cellfun(@(s) sprintf('''%s''', s), fullTableNames, 'uni', false);
    inList = strjoin(quoted, ',');

    sql = [ ...
        'SELECT CONCAT(''`'',TABLE_SCHEMA,''`.'',''`'',TABLE_NAME,''`'') AS child, ' ...
        'CONCAT(''`'',REFERENCED_TABLE_SCHEMA,''`.'',''`'',REFERENCED_TABLE_NAME,''`'') AS parent ' ...
        'FROM INFORMATION_SCHEMA.KEY_COLUMN_USAGE ' ...
        'WHERE REFERENCED_TABLE_SCHEMA IS NOT NULL ' ...
        'AND CONCAT(''`'',REFERENCED_TABLE_SCHEMA,''`.'',''`'',REFERENCED_TABLE_NAME,''`'') IN (' inList ')'];

    try
        res = dj.conn().query(sql);
    catch
        return; % If the metadata query fails, fall back to no external tables.
    end

    if isempty(res)
        return;
    end

    childTables = res.child;
    for idx = 1:numel(childTables)
        clsName = dj.conn().tableToClass(childTables{idx}, false);
        try
            cls = feval(clsName);
            extList = [extList, cls.descendants]; %#ok<AGROW>
        catch
            % If the class cannot be resolved/loaded, skip it quietly.
        end
    end

    extList = unique(extList, 'stable');
end