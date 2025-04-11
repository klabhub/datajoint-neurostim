function nsParpopulate(tbl,cls,opts,pv)
% Run parpopulate on as many workers as there are rows to fill.
%
% EXAMPLE
% opts =  {'time','00:30:00','mem','16GB','cpus-per-task',1,'requeue',''}
% restrict = {'p<0.05','subject="joe"'} % define restrictions as a cell
%           array of chars. 
% nsParPopulate(cr.Crf,cls,opts,restrict=restrict);
arguments
    tbl (1,1) dj.Relvar     % table to populate
    cls (1,1) mslurm        % slurm cluster to use
    opts (1,:) cell = {'time','00:30:00','mem','16GB','cpus-per-task',1,'requeue',''}; % sbatch options for mslurm
    pv.maxWorkers = 64    % Never more than this number of workers
    pv.restrict (1,:) cell = {}  % Restrict parpopulate(tbl, restrict{:})  
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
nrToDo = count(getKeySource(tbl) - tbl);
nrWorkers = min(nrToDo,pv.maxWorkers);
if nrToDo >0
    fprintf("Populating %d rows of %s with %d workers (cmd=%s)\n",nrToDo,tbl.className,nrWorkers,cmd);
    cls.remote(cmd,'nrWorkers',nrWorkers,'expressionName',exp,'sbatchOptions',opts); 
else
    fprintf("Nothing to populate for %s (table has %d rows)\n",cmd,count(tbl));
end

