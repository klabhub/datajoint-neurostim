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

if isa(tbl,'table')
    % User passed the output of jobs, which is a matlab table
    tpl =table2struct(tbl(row,:));
else
    tpl = fetch(tbl ,'*',sprintf('LIMIT 1 OFFSET %d',row-1));
end

%% Setup the break point
if isempty(tpl.error_stack)
    fprintf('No error stack; just retrying without a breakpoint\n');
else
    % Determine the function that caused the error
    ffile = tpl.error_stack(1).file;
    if contains(ffile,'+')
        % Part of a (nested) namespace - resolve to the a.b.c name
        if contains(ffile,'/')
            fs = '/';
        else
            fs = '\';
        end
        fname =extractAfter(ffile,[fs '+']);
        fname  = strrep(fname,'+','');
        fname =extractBefore(fname,'.m');
        fname = strrep(fname,fs,'.');
    else
        % No namespace used: fname should be on the path.
        [~,fname] =fileparts(ffile);
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
if ~isa(tbl,'table')
    T = ns.Jobs & struct('table_name',tpl.table_name,'error_message',tpl.error_message);
end