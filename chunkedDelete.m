function chunkedDelete(targetQuery, batchSize)
    if nargin < 2; batchSize = 100; end

    % 1. Dry-Run Analysis
    disp('Analyzing dependency tree and tuple counts...');
    list = targetQuery.descendants;
    list = flip(list); 
    
    totalTuples = 0;
    summary = cell(length(list), 2);
    for i = 1:length(list)
        tableObj = feval(list{i});
        toDelete = tableObj & targetQuery.restriction;
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
        toDelete = tableObj & targetQuery.restriction;
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