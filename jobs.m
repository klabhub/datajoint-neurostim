function keyT = jobs (tbl,pv)
% Utility funciton to get information on the information in the jobs tabel
% Mainly intended to quickly show which file (errfile) generated which
% error (message) for which item.
%  INPUT
% 'status'   =Cell string of jobs with a specific status to include {'error'}
arguments
    tbl (1,1) = ns.Jobs
    pv.status = {'error'}
end
if count(tbl)==0;fprintf('No jobs in the %s table\n',tbl.className);return;end
[keys,tablename,status, message,stack,timestamp] = fetchn(tbl &struct('status',pv.status),'key','table_name','status','error_message','error_stack','timestamp');
keyT  = struct2table(catstruct(1,keys{:})); % Add all elements of the key as columns to the table
% Extract the top of the error stack (when relevant)
errfile = repmat("",[height(keyT) 1]);
for i=1:numel(stack)
    if isstruct(stack{i})
        errfile(i) = string(stack{i}(1).file);
    end
end
% Combine and sort
keyT =addvars(keyT,tablename,status,message,errfile, stack,timestamp);
keyT = convertvars(keyT,@iscellstr,'string');
keyT = movevars(keyT,["errfile","message"],"Before",1);
keyT = sortrows(keyT,"timestamp","descend");
end