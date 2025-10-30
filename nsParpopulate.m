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
%
% Use dryrun = true to get a brief report of the jobs that would be started
% with dryrun =false.
%
arguments
    tbl (1,1) dj.Relvar     % table to populate
    cls (1,1) mslurm        % slurm cluster to use
    opts (1,:) cell = {'time','00:30:00','mem','16GB','cpus-per-task',1,'requeue',''}; % sbatch options for mslurm
    pv.maxWorkers = 64    % Never more than this number of workers
    pv.restrict (1,:)  = {}  % Restrict parpopulate(tbl, restrict{:})
    pv.clearJobStatus (1,:) string = "error" %  string array of status flags that should be cleared from the jobs table.
    pv.dryrun (1,1) logical = false % Set to true to get command line feedback on which jobs would be started
    pv.env (1,:) string = "" % Additions to the environment on the cluster for this call only
    pv.jobName (1,1) string = "" % Informative name for the job used in SLURM sacct. Does not need to be unique
end
warnState = warning('query');
warning('off','DataJoint:longCondition');
if pv.dryrun
    drMsg = "[DRYRUN]";
else
    drMsg = "";
end
if ~isempty(tbl.restrictions)
    fprintf(2,' The restriction currently specified on the table will be ignored. (You probably want to use pv.restrict )\n');
end
if ischar(pv.restrict) || isstring(pv.restrict)
    pv.restrict= {pv.restrict};
end
if numel(pv.restrict)>0
    % Restrictions - each must be a string
    assert(all(cellfun(@ischar,pv.restrict)),'Each element of the restriction must be a char')
    restriction = ['''' strjoin(pv.restrict,''',''') ''''];
    cmd = sprintf("parpopulate(%s,%s)",tbl.className,restriction);
else
    % Unrestricted populate
    cmd = sprintf("parpopulate(%s)",tbl.className);
end
% Construct a name for this expression (has to be a valid script name in
% matlab)
exp = sprintf("parpop_%s",strrep(tbl.className,'.','_'));

keysource =(tbl.getKeySource & pv.restrict) -tbl;

if exists(keysource) &&  pv.clearJobStatus~=""
    % Check the jobs table to see if something needs to be cleared
    jobsTableName = extractBefore(tbl.className,'.') + ".Jobs";
    if exist(jobsTableName,"class")
        % Find related jobs (i.e. filling the same table and the specified status)
        sameTableJobs = feval(jobsTableName) & in("status",pv.clearJobStatus) & sprintf('table_name="%s"',tbl.className);
        if count(sameTableJobs)>0
            % Determine which of these jobs will be retried by the current
            % populate call
            keysourceTable = fetchtable(keysource);
            % Get the complete jobs table as a matlab table and join with
            % the keysrouce to get the jobs that will be retried
            retriedJobs = innerjoin(jobs(sameTableJobs),keysourceTable,'Keys',keysourceTable.Properties.VariableNames,'LeftVariables',["table_name" "key_hash"]);
            fprintf("%s Deleting %d jobs from %s to retry\n",drMsg,height(retriedJobs),jobsTableName)
            if ~pv.dryrun
                delQuick( sameTableJobs & table2struct(retriedJobs) )
            end
        end
    else
        fprintf('No jobs table (%s) on disk. \n',jobsTableName)
    end
end



nrChildren = count(tbl);
nrToDo = count(keysource);
nrWorkers = min(nrToDo,pv.maxWorkers);
if nrToDo >0
    fprintf("%s Populating %d rows of %s with %d workers (cmd=%s)\n",drMsg,nrToDo,tbl.className,nrWorkers,cmd);
    if ~pv.dryrun
        if pv.env ==""
            env = pv.env;
        else
            env = [cls.env pv.env];
        end
        cls.remote(cmd,'nrWorkers',nrWorkers,'expressionName',exp,'sbatchOptions',opts,'env',env,'jobName',pv.jobName);
    end
else
    fprintf("%s Nothing to populate for %s (table has %d rows already computed)\n",drMsg,cmd,nrChildren);
end

warning(warnState);