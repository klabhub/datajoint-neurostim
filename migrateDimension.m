function migrateDimension(dbName, startAt)
% migrateDimension  Migrate ns.Dimension and dependents to the new
%   design with a paradigm-specific primary key and separate DimensionParm table.
% This design was introduced in May 2026 to address the problem that the old way of using 
% define() did not store all of the parameters needed to reconstruct the dimension conditions, 
% and also did not allow for dimensions to be populated automatically upon insertion of anew experiments.
%
%   migrateDimension(dbName) runs all four phases in sequence:
%     1. backup   - Copy Tuning, Epoch, EpochChannel to _bak tables,
%                   adding 'paradigm' from ns/experiment.
%     2. drop     - Drop ns.Dimension and all dependent tables (after
%                   confirmation).
%     3. recreate - Recreate tables via DataJoint schema definitions and
%                   verify 'paradigm' is in the new ns.Dimension PK.
%     4. restore  - INSERT Epoch, EpochChannel, and Tuning data from _bak.
%
%   migrateDimension(dbName, startAt) starts from the specified phase,
%   skipping earlier phases.  Use this to resume after an interruption.
%
%   NOTE: ns/#dimension_parm must already exist with PK (dimension, paradigm)
%   before running this function.
%
%   EXAMPLES
%     migrateDimension('djEEG')             % run all phases
%     migrateDimension('djEEG', 'drop')     % resume from drop phase
%     migrateDimension('djEEG', 'restore')  % only restore from backups

arguments
    dbName  (1,1) string
    startAt (1,1) string ...
        {mustBeMember(startAt, {'backup','drop','recreate','restore'})} = 'backup'
end

phases   = {'backup','drop','recreate','restore'};
startIdx = find(strcmpi(startAt, phases));
fprintf('=== migrateDimension: starting from phase ''%s'' (database: %s) ===\n\n', ...
    startAt, dbName);

conn = dj.conn();
db   = char(dbName);

% ==================================================================
% PRE-FLIGHT: verify ns/#dimension_parm and coverage
% ==================================================================
checkDimensionParm(conn, db);

% ==================================================================
% PHASE 1: BACKUP
% ==================================================================
if startIdx <= 1
    fprintf('--- Phase 1: Backup ---\n');
    tblTuning    = findTable(conn, db, '%tuning',        'NOT LIKE ''%bak%''',           'Tuning');
    tblEpoch     = findTable(conn, db, '%epoch',         'NOT LIKE ''%bak%'' AND table_name NOT LIKE ''%channel%''', 'Epoch');
    tblEpochChan = findTable(conn, db, '%epoch%channel', 'NOT LIKE ''%bak%''',           'EpochChannel');
    backupWithParadigm(conn, db, tblTuning,    [tblTuning    '_bak']);
    backupWithParadigm(conn, db, tblEpoch,     [tblEpoch     '_bak']);
    backupWithParadigm(conn, db, tblEpochChan, [tblEpochChan '_bak']);
    fprintf('Backups complete.\n\n');
end

% ==================================================================
% PHASE 2: DROP
% ==================================================================
if startIdx <= 2
    fprintf('--- Phase 2: Drop ---\n');
    tblTuning    = findTable(conn, db, '%tuning',        'NOT LIKE ''%bak%''',           'Tuning');
    tblEpoch     = findTable(conn, db, '%epoch',         'NOT LIKE ''%bak%'' AND table_name NOT LIKE ''%channel%''', 'Epoch');
    tblEpochChan = findTable(conn, db, '%epoch%channel', 'NOT LIKE ''%bak%''',           'EpochChannel');
    tblDimension = findTableOptional(conn, db, 'ns/dimension',          'NOT LIKE ''%bak%'' AND table_name NOT LIKE ''%condition%'' AND table_name NOT LIKE ''%trial%''', 'Dimension');
    tblDimCond   = findTableOptional(conn, db, '%dimension%condition%', 'NOT LIKE ''%bak%''', 'DimensionCondition');
    tblDimTrial  = findTableOptional(conn, db, '%dimension%trial%',     'NOT LIKE ''%bak%''', 'DimensionTrial');
    fprintf('The following tables will be DROPPED:\n');
    fprintf('  %s\n', tblEpochChan);
    fprintf('  %s\n', tblEpoch);
    fprintf('  %s\n', tblTuning);
    if ~isempty(tblDimCond),  fprintf('  %s\n', tblDimCond);  end
    if ~isempty(tblDimTrial), fprintf('  %s\n', tblDimTrial); end
    if ~isempty(tblDimension), fprintf('  %s\n', tblDimension); end

    answer = input('\nType ''yes'' to proceed with DROP: ', 's');
    if ~strcmpi(strtrim(answer), 'yes')
        fprintf('Aborted. No tables were dropped.\n');
        return;
    end

    conn.query('SET FOREIGN_KEY_CHECKS = 0');
    dropTable(conn, db, tblEpochChan);
    dropTable(conn, db, tblEpoch);
    dropTable(conn, db, tblTuning);
    if ~isempty(tblDimCond),   dropTable(conn, db, tblDimCond);   end
    if ~isempty(tblDimTrial),  dropTable(conn, db, tblDimTrial);  end
    if ~isempty(tblDimension), dropTable(conn, db, tblDimension); end
    conn.query('SET FOREIGN_KEY_CHECKS = 1');
    fprintf('All tables dropped.\n\n');
end

% ==================================================================
% PHASE 3: RECREATE
% ==================================================================
if startIdx <= 3
    fprintf('--- Phase 3: Recreate via DataJoint ---\n');
    % Flush DataJoint's schema cache so it re-reads table definitions
    % from MySQL rather than using the state from before the drops.
    ns.Dimension.schema.reload(true);
    createViaDataJoint('ns.Dimension',          @() count(ns.Dimension));
    createViaDataJoint('ns.DimensionCondition', @() count(ns.DimensionCondition));
    createViaDataJoint('ns.DimensionTrial',     @() count(ns.DimensionTrial));
    createViaDataJoint('ns.Epoch',              @() count(ns.Epoch));
    createViaDataJoint('ns.EpochChannel',       @() count(ns.EpochChannel));
    createViaDataJoint('ns.Tuning',             @() count(ns.Tuning));

    % Re-discover dimension table name now that it has been recreated
    tblDimension = findTable(conn, db, 'ns/%dimension', ...
        'NOT LIKE ''%bak%'' AND table_name NOT LIKE ''%condition%'' AND table_name NOT LIKE ''%trial%''', ...
        'Dimension');

    fprintf('\nVerifying primary key of ns.Dimension includes ''paradigm''...\n');
    pkCols = conn.query(sprintf([ ...
        'SELECT column_name FROM information_schema.key_column_usage ' ...
        'WHERE table_schema = ''%s'' AND table_name = ''%s'' ' ...
        'AND constraint_name = ''PRIMARY'' ' ...
        'ORDER BY ordinal_position'], db, tblDimension));
    fprintf('  Primary key: %s\n', strjoin(pkCols.column_name, ', '));
    if any(strcmpi(pkCols.column_name, 'paradigm'))
        fprintf('  OK: ''paradigm'' is in the primary key.\n\n');
    else
        warning('migrateDimension:noPKParadigm', ...
            '''paradigm'' is NOT in the primary key. Check ns.Dimension schema.');
    end

    fprintf('Populating ns.Dimension...\n');
    parpopulate(ns.Dimension);
    fprintf('ns.Dimension populate complete (%d rows).\n\n', count(ns.Dimension));
end

% ==================================================================
% PHASE 4: RESTORE
% ==================================================================
fprintf('--- Phase 4: Restore from backups ---\n');
tblTuning    = findTable(conn, db, '%tuning',        'NOT LIKE ''%bak%''',           'Tuning');
tblEpoch     = findTable(conn, db, '%epoch',         'NOT LIKE ''%bak%'' AND table_name NOT LIKE ''%channel%''', 'Epoch');
tblEpochChan = findTable(conn, db, '%epoch%channel', 'NOT LIKE ''%bak%''',           'EpochChannel');
tblDimension = findTable(conn, db, '%dimension', ...
    'NOT LIKE ''%bak%'' AND table_name NOT LIKE ''%condition%'' AND table_name NOT LIKE ''%trial%''', ...
    'Dimension');
restoreFromBackup(conn, db, [tblEpoch     '_bak'], tblEpoch,     tblDimension);
restoreFromBackup(conn, db, [tblEpochChan '_bak'], tblEpochChan, tblDimension);
restoreFromBackup(conn, db, [tblTuning    '_bak'], tblTuning,    tblDimension);
fprintf('Restore complete.\n\n');

fprintf('=== migrateDimension complete (database: %s) ===\n', dbName);
end


% ======================================================================
function checkDimensionParm(conn, db)
% Verify ns/#dimension_parm exists, then find any (paradigm, dimension)
% pairs present in ns/dimension but missing from ns/#dimension_parm.

% 1. Check that the DimensionParm table exists
parmTbl = conn.query(sprintf([ ...
    'SELECT table_name FROM information_schema.tables ' ...
    'WHERE table_schema = ''%s'' AND table_name LIKE ''%%dimension_parm%%'' ' ...
    'AND table_name NOT LIKE ''%%paradigm%%'' LIMIT 1'], db));
if isempty(parmTbl.table_name)
    error('migrateDimension:noDimensionParm', ...
        'Table ns/#dimension_parm not found in ''%s''. ', ...
        'Create it with the new (dimension, paradigm) PK before migrating.', db);
end

% 2. Check that ns/dimension exists (may not exist if starting late)
dimTbl = conn.query(sprintf([ ...
    'SELECT table_name FROM information_schema.tables ' ...
    'WHERE table_schema = ''%s'' AND table_name = ''ns/dimension'''], db));
if isempty(dimTbl.table_name)
    fprintf('Pre-flight: ns/dimension not found — skipping coverage check.\n\n');
    return
end

% Check that all Dimensions are now defined in DimensionParm
dimsWithoutParms  = (ns.Dimension * proj(ns.Experiment,'paradigm')) - proj(ns.DimensionParm);

if exists(dimsWithoutParms)
    fprintf('Thse dimensions are missing from DimensionParm\n')
    disp(dimsWithoutParms)
else 
    fprintf('All dimensions have a parm in DimensionParm\n')

end
end

% ======================================================================
function name = findTable(conn, db, pattern, excludeClause, label)
% Find a table matching LIKE pattern, with an optional exclusion clause.
result = conn.query(sprintf([ ...
    'SELECT table_name FROM information_schema.tables ' ...
    'WHERE table_schema = ''%s'' AND table_name LIKE ''%s'' ' ...
    'AND table_name %s ' ...
    'ORDER BY table_name LIMIT 1'], db, pattern, excludeClause));
if isempty(result) || isempty(result.table_name)
    error('migrateDimension:tableNotFound', ...
        'No table matching ''%s'' (%s) found in ''%s''.', pattern, label, db);
end
name = result.table_name{1};
end


% ======================================================================
function name = findTableOptional(conn, db, pattern, excludeClause, label)
% Like findTable but returns '' instead of erroring if not found.
try
    name = findTable(conn, db, pattern, excludeClause, label);
catch
    fprintf('  (optional table not found, skipping: %s)\n', label);
    name = '';
end
end


% ======================================================================
function backupWithParadigm(conn, db, srcTable, dstTable)
existing = conn.query(sprintf([ ...
    'SELECT table_name FROM information_schema.tables ' ...
    'WHERE table_schema = ''%s'' AND table_name = ''%s'''], db, dstTable));
if ~isempty(existing.table_name)
    fprintf('  Backup already exists, skipping: %s\n', dstTable);
    return;
end
fprintf('  Creating backup: %s ...', dstTable);
conn.query(sprintf([ ...
    'CREATE TABLE `%s`.`%s` AS ' ...
    'SELECT t.*, e.paradigm ' ...
    'FROM `%s`.`%s` t ' ...
    'JOIN `%s`.`ns/experiment` e ' ...
    '  ON t.starttime    = e.starttime ' ...
    ' AND t.session_date = e.session_date ' ...
    ' AND t.subject      = e.subject'], ...
    db, dstTable, db, srcTable, db));
cnt    = conn.query(sprintf('SELECT COUNT(*) AS n FROM `%s`.`%s`', db, dstTable));
cntSrc = conn.query(sprintf('SELECT COUNT(*) AS n FROM `%s`.`%s`', db, srcTable));
fprintf(' %d / %d rows', cnt.n, cntSrc.n);
if cnt.n == cntSrc.n
    fprintf(' - OK\n');
else
    fprintf(' - WARNING: count mismatch!\n');
end
end


% ======================================================================
function dropTable(conn, db, tblName)
conn.query(sprintf('DROP TABLE IF EXISTS `%s`.`%s`', db, tblName));
fprintf('  Dropped: %s\n', tblName);
end


% ======================================================================
function restoreFromBackup(conn, db, srcTable, dstTable, dimTable)
fprintf('  Restoring %s ...', dstTable);
dstCols = conn.query(sprintf([ ...
    'SELECT column_name FROM information_schema.columns ' ...
    'WHERE table_schema = ''%s'' AND table_name = ''%s'' ' ...
    'ORDER BY ordinal_position'], db, dstTable));
srcCols = conn.query(sprintf([ ...
    'SELECT column_name FROM information_schema.columns ' ...
    'WHERE table_schema = ''%s'' AND table_name = ''%s'' ' ...
    'ORDER BY ordinal_position'], db, srcTable));
commonCols = dstCols.column_name( ...
    ismember(lower(dstCols.column_name), lower(srcCols.column_name)));
if isempty(commonCols)
    error('migrateDimension:noCommonColumns', ...
        'No common columns between %s and %s.', srcTable, dstTable);
end
insertCols = strjoin(cellfun(@(c) sprintf('`%s`', c), commonCols, ...
    'UniformOutput', false), ', ');
selectCols = strjoin(cellfun(@(c) sprintf('s.`%s`', c), commonCols, ...
    'UniformOutput', false), ', ');

% Build JOIN condition on the PK of the dimension table so that only rows
% with a matching entry in the new Dimension table are copied.
pkRes = conn.query(sprintf([ ...
    'SELECT column_name FROM information_schema.key_column_usage ' ...
    'WHERE table_schema = ''%s'' AND table_name = ''%s'' ' ...
    'AND constraint_name = ''PRIMARY'' ' ...
    'ORDER BY ordinal_position'], db, dimTable));
joinCond = strjoin(cellfun(@(c) sprintf('s.`%s` = d.`%s`', c, c), ...
    pkRes.column_name, 'UniformOutput', false), ' AND ');

conn.query(sprintf([ ...
    'INSERT IGNORE INTO `%s`.`%s` (%s) ' ...
    'SELECT %s FROM `%s`.`%s` s ' ...
    'JOIN `%s`.`%s` d ON %s'], ...
    db, dstTable, insertCols, selectCols, db, srcTable, db, dimTable, joinCond));
cntSrc = conn.query(sprintf('SELECT COUNT(*) AS n FROM `%s`.`%s`', db, srcTable));
cnt    = conn.query(sprintf('SELECT COUNT(*) AS n FROM `%s`.`%s`', db, dstTable));
fprintf(' %d / %d rows done.\n', cnt.n, cntSrc.n);
end


% ======================================================================
function createViaDataJoint(label, accessFn)
fprintf('  %-30s ... ', label);
try
    accessFn();
    fprintf('OK\n');
catch ME
    fprintf('ERROR: %s\n', ME.message);
end
end
