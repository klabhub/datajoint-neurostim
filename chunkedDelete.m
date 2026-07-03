function chunkedDelete(targetQuery, batchSize, extraRelvars, autoDetectFK)
% To avoid long server locks, delete rows in tables in batches.
% This will first scan the query tree, and then delete rows from the bottom
% up to avoid long table locks.
% Batch size is adjusted per table from estimated row size when metadata is
% available, and falls back to the provided row-count batchSize otherwise.
% If dj.congif('safeMode') =1 the function asks for confirmation (once, after showing
% the deletion plan)/
arguments
    targetQuery (1,1) dj.internal.GeneralRelvar
    batchSize (1,1) double = 500 %Delete this many rows each call.
    extraRelvars cell = {} % Optional cell array of additional relvars (possibly cross-schema) to include.
    autoDetectFK (1,1) logical = true % When true, include external FK children discovered via INFORMATION_SCHEMA
end
wrnStatus = warning("query");
warning("off",'DataJoint:longCondition')
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

referenceRowBytes = NaN;
if ~isempty(descBase)
    referenceRowBytes = estimateRowBytes(feval(descBase{1}).fullTableName);
end

list = flip(list);

targetBytesPerDelete = resolveDeleteByteTarget(batchSize, referenceRowBytes);
maxUnchunkedBytes = resolveMaxUnchunkedDeleteBytes();


totalTuples = 0;
summary = cell(length(list), 6);
for i = 1:length(list)
    tableObj = feval(list{i});
    toDelete = tableObj & targetQuery.proj();
    c = count(toDelete);
    [rowBytes, adaptiveBatchSize] = estimateDeleteBatchSize(tableObj, batchSize, targetBytesPerDelete);
    isLeafTable = isLeafDeleteTable(tableObj, autoDetectFK);
    useSingleShotDelete = isLeafTable && shouldDeleteInSingleShot(c, rowBytes, maxUnchunkedBytes);
    summary{i, 1} = list{i};
    summary{i, 2} = c;
    summary{i, 3} = rowBytes;
    summary{i, 4} = adaptiveBatchSize;
    summary{i, 5} = isLeafTable;
    summary{i, 6} = useSingleShotDelete;
    totalTuples = totalTuples + c;
end


% Display Plan
fprintf('\n--- DELETE PLAN (Bottom-Up) ---\n');
for i = 1:size(summary, 1)
    if summary{i, 2} > 0
        modeLabel = 'chunked';
        if summary{i, 6}
            modeLabel = 'single-shot leaf';
        end
        if isnan(summary{i, 3})
            fprintf('%8d tuples from %s (%s, batch %d, row size unavailable)\n', ...
                summary{i, 2}, summary{i, 1}, modeLabel, summary{i, 4});
        else
            fprintf('%8d tuples from %s (%s, batch %d, %.1f KB/row)\n', ...
                summary{i, 2}, summary{i, 1}, modeLabel, summary{i, 4}, summary{i, 3}/1024);
        end
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

    adaptiveBatchSize = summary{i, 4};
    useSingleShotDelete = summary{i, 6};

    if useSingleShotDelete
        fprintf('\nProcessing %s with single-shot delete...\n', tableName);
    else
        fprintf('\nProcessing %s with batch size %d...\n', tableName, adaptiveBatchSize);
    end
    keys = fetch(toDelete); % Fetch keys into MATLAB RAM

    if useSingleShotDelete
        delQuick(tableObj & keys);
        fprintf('  Deleted %d/%d\n', numTuples, numTuples);
        continue;
    end

    j = 1;
    while j <= numTuples
        endIdx = min(j + adaptiveBatchSize - 1, numTuples);
        currentBatch = keys(j:endIdx);

        try
            delQuick(tableObj & currentBatch); % Commit-per-batch
            fprintf('  Deleted %d/%d\n', endIdx, numTuples);
            j = j + adaptiveBatchSize; % Advance only on success
        catch ME
            if contains(ME.message, 'gone away') || contains(ME.message, 'Lost connection')
                fprintf('\n[CONNECTION LOST] Attempting to reconnect...\n');
                pause(5); % Give the server a moment to finish restarting
                try
                    dj.conn().reopen(); % Re-establish the pipe
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
warning(wrnStatus)
end


function bytesPerDelete = resolveDeleteByteTarget(fallbackBatchSize, referenceRowBytes)
% Keep a roughly constant amount of table data per delete batch.

    bytesPerDelete = getenv("NS_BYTESPERDELETE");
    if isempty(bytesPerDelete)
        if ~isnan(referenceRowBytes) && referenceRowBytes > 0
            bytesPerDelete = fallbackBatchSize * referenceRowBytes;
        else
            bytesPerDelete = fallbackBatchSize * 8192;
        end
        return;
    end

    bytesPerDelete = str2double(bytesPerDelete);
    if isnan(bytesPerDelete) || bytesPerDelete <= 0
        error('NS_BYTESPERDELETE must be a positive numeric value.');
    end
end


function maxUnchunkedBytes = resolveMaxUnchunkedDeleteBytes()
% Limit when a leaf table may be deleted in a single statement.

    maxUnchunkedBytes = getenv("NS_MAXUNCHUNKEDDELETE");
    if isempty(maxUnchunkedBytes)
        maxUnchunkedBytes = 50e6;
        return;
    end

    maxUnchunkedBytes = str2double(maxUnchunkedBytes);
    if isnan(maxUnchunkedBytes) || maxUnchunkedBytes <= 0
        error('NS_MAXUNCHUNKEDDELETE must be a positive numeric value.');
    end
end


function [rowBytes, adaptiveBatchSize] = estimateDeleteBatchSize(tableObj, fallbackBatchSize, targetBytesPerDelete)
% Use INFORMATION_SCHEMA estimates when possible, otherwise keep the input batch size.

    rowBytes = estimateRowBytes(tableObj.fullTableName);
    if isnan(rowBytes) || rowBytes <= 0
        adaptiveBatchSize = fallbackBatchSize;
        return;
    end

    adaptiveBatchSize = max(1, floor(targetBytesPerDelete / rowBytes));
end


function tf = shouldDeleteInSingleShot(numTuples, rowBytes, maxUnchunkedBytes)
% Only unchunk leaf deletes when the estimated payload is small.

    tf = ~isnan(rowBytes) && rowBytes > 0 && (numTuples * rowBytes) <= maxUnchunkedBytes;
end


function tf = isLeafDeleteTable(tableObj, includeExternalChildren)
% A leaf delete table has no dependent child tables.

    directChildren = normalizeTableNameList(tableObj.children());
    if includeExternalChildren
        externalChildren = normalizeTableNameList(discoverDirectExternalChildren(tableObj.fullTableName));
        directChildren = unique([directChildren, externalChildren], 'stable');
    end
    tf = isempty(directChildren);
end


function names = normalizeTableNameList(names)
% Convert table-name outputs to a row cell array for safe concatenation.

    if isempty(names)
        names = {};
    elseif ischar(names)
        names = {names};
    elseif isstring(names)
        names = cellstr(names);
    end

    names = reshape(names, 1, []);
end


function rowBytes = estimateRowBytes(fullTableName)
% Estimate bytes per stored row using table metadata.

    rowBytes = NaN;
    [schemaName, tableName] = splitFullTableName(fullTableName);
    sql = sprintf([ ...
        'SELECT AVG_ROW_LENGTH, DATA_LENGTH, TABLE_ROWS ' ...
        'FROM INFORMATION_SCHEMA.TABLES ' ...
        'WHERE TABLE_SCHEMA = ''%s'' AND TABLE_NAME = ''%s'''], ...
        schemaName, tableName);

    try
        res = dj.conn().query(sql);
    catch
        return;
    end

    if isempty(res)
        return;
    end

    if isfield(res, 'AVG_ROW_LENGTH') && ~isempty(res.AVG_ROW_LENGTH) && ~isnan(res.AVG_ROW_LENGTH)
        rowBytes = double(res.AVG_ROW_LENGTH);
        if rowBytes > 0
            return;
        end
    end

    if isfield(res, 'DATA_LENGTH') && isfield(res, 'TABLE_ROWS') && ...
            ~isempty(res.TABLE_ROWS) && res.TABLE_ROWS > 0 && ...
            ~isempty(res.DATA_LENGTH) && ~isnan(res.DATA_LENGTH)
        rowBytes = double(res.DATA_LENGTH) / double(res.TABLE_ROWS);
    end
end


function [schemaName, tableName] = splitFullTableName(fullTableName)
% Convert `schema`.`table` into unquoted names for INFORMATION_SCHEMA.

    parts = regexp(fullTableName, '^`([^`]+)`\.`([^`]+)`$', 'tokens', 'once');
    if isempty(parts)
        error('Could not parse table name %s.', fullTableName);
    end

    schemaName = parts{1};
    tableName = parts{2};
end


function childTables = discoverDirectExternalChildren(fullTableName)
% Query INFORMATION_SCHEMA for direct FK children across all schemas.

    childTables = {};
    [schemaName, tableName] = splitFullTableName(fullTableName);
    sql = sprintf([ ...
        'SELECT CONCAT(''`'',TABLE_SCHEMA,''`.'',''`'',TABLE_NAME,''`'') AS child ' ...
        'FROM INFORMATION_SCHEMA.KEY_COLUMN_USAGE ' ...
        'WHERE REFERENCED_TABLE_SCHEMA = ''%s'' ' ...
        'AND REFERENCED_TABLE_NAME = ''%s'''], ...
        schemaName, tableName);

    try
        res = dj.conn().query(sql);
    catch
        return;
    end

    if isempty(res) || ~isfield(res, 'child')
        return;
    end

    childTables = unique(res.child, 'stable');
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