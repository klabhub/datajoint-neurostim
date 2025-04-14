function nsParpopulate(tbl,cls,opts,pv)
% Run parpopulate on as many workers as there are rows to fill.
%
% EXAMPLE
% opts =  {'time','00:30:00','mem','16GB','cpus-per-task',1,'requeue',''}
% restrict = {'p<0.05','subject="joe"'} % define restrictions as a cell
%           array of chars.
% nsParPopulate(cr.Crf,cls,opts,restrict=restrict);
%
% By default, error jobs in the .Jobs table will be deleted before posting
% the new jobs. You can disable this with clearJobStatus = "", or you can
% also delete reserved jobs with clearJobStatus = ["error" "reserved"]
% (if those jobs are currently running this will cause problems...)
%
% The maximum number of workers is set to 64 by default. Be nice.
arguments
    tbl (1,1) dj.Relvar     % table to populate
    cls (1,1) mslurm        % slurm cluster to use
    opts (1,:) cell = {'time','00:30:00','mem','16GB','cpus-per-task',1,'requeue',''}; % sbatch options for mslurm
    pv.maxWorkers = 64    % Never more than this number of workers
    pv.restrict (1,:) cell = {}  % Restrict parpopulate(tbl, restrict{:})
    pv.clearJobStatus (1,:) string = "error" %  string array of status flags that should be cleared from the jobs table.
    pv.dryrun (1,1) logical = false % Set to true to get command line feedback on which jobs would be started
end
if pv.dryrun 
    drMsg = "[DRYRUN]";
else
    drMsg = "";
end

if numel(pv.restrict)>0
    % Restrictions - each must be a string
    assert(all(cellfun(@ischar,pv.restrict)),'Each element of the restriction must be a char')
    r = ['''' strjoin(pv.restrict,''',''') ''''];
    cmd = sprintf("parpopulate(%s,%s)",tbl.className,r);
else
    % Unrestricted populate
    cmd = sprintf("parpopulate(%s)",tbl.className);
end
% Construct a name for this expression (has to be a valid script name in
% matlab)
exp = sprintf("parpop_%s",strrep(tbl.className,'.','_'));

if pv.clearJobStatus~=""
    % Check the jobs table to see if something needs to be cleared
    jobsTable = extractBefore(tbl.className,'.') + ".Jobs";
    jt = feval(jobsTable);
    jt = jt & in("status",pv.clearJobStatus) & sprintf('table_name="%s"',tbl.className);
    if count(jt)>0
        fprintf('%s %d related jobs in %s will be deleted',drMsg,count(jt),jobsTable)
        if ~pv.dryrun
            delQuick(jt)
        end
    end
end

nrToDo = count(getKeySource(tbl) - proj(tbl));
nrWorkers = min(nrToDo,pv.maxWorkers);
if nrToDo >0
    fprintf("%s Populating %d rows of %s with %d workers (cmd=%s)\n",drMsg,nrToDo,tbl.className,nrWorkers,cmd);
    if ~pv.dryrun
        cls.remote(cmd,'nrWorkers',nrWorkers,'expressionName',exp,'sbatchOptions',opts);
    end
else
    fprintf("%s Nothing to populate for %s (table has %d rows)\n",drMsg,cmd,count(tbl));
end

