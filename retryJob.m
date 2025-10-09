function T= retryJob(row,tbl)
% Retry populating an item in the Jobs table that caused an error. This
% retry will put a breakpoint on the line where the previous attempt
% failed and stop execution there. 
% 
% If this funciton runs to completion, and is successful, the job will be
% removed from the jobs table.
% 
% Returns the subset of the Jobs table with the same error message (and
% table_name). If the problem is fixed, the user can call del on this
% returen table to clear those errors.
arguments
    row  = 1  % The row of the Jobs table to process
    tbl = ns.Jobs
end

tpl = fetch(tbl ,'*',sprintf('LIMIT 1 OFFSET %d',row-1));

%% Setup the break point
if isempty(tpl.error_stack)
    fprintf('No error stack; just retrying without a breakpoint\n');
else
[p,f] =fileparts(tpl.error_stack(1).file);
[~,nmspc] =fileparts(p);
if contains(nmspc,'+')
    fname = [nmspc(2:end) '.' f];
else
    fname =f;
end
fprintf('Setting breakpoint at line %d in %s\n',tpl.error_stack(1).line,fname);
dbstop('in', string(fname),'at', string(tpl.error_stack(1).line));
end

%% Call populate
try
    populate(feval(tpl.table_name),tpl.key);
    success=true;
catch me
    me.message
    success=false;    
end

%% Success - remove the item from the jobs table
if success
    delQuick(ns.Jobs & ns.stripToPrimary(ns.Jobs,tpl))
end

T = ns.Jobs & struct('table_name',tpl.table_name,'error_message',tpl.error_message);