function pool = nsParPool(pv)
% Function to start (or not start) a local parpool based on an 
% environment variable. 
% 
% This is useful to avoid storing or hardcoding things like nrWorkers 
% in a table for a makeTuples function  (which forces the same number 
% of workers forever and is annoying when debugging).
% 
% To use this, write parfor code like this
% 
% pool = nsParpool;
% if isempty(pool)
%   for i=1:N
%       doIt(i);
%   end
% else 
%   parfor (i=1:N)
%       doIt(i);
%   end
% end
%
% If you want to use the parfor with 10 workers, call setenv("NS_PARPOOL",10) 
% or you start a parpool directly with a call to parpool() and ensure that
% NS_PARPOOL environment variable is not set (or "")
%
% To temporarily use the for-loop even with a parpool running, call
% setenv("NS_PARPOOL",0)
%
% If this function creates the pool, SpmdEnabled is set to false by
% default. Set the spmd parameter to true to enable it on the pool.
%
% Calling this function with the nrWorkers argument set bypassess the
% environment variable (you can use this to start a pool manually for
% instance, but using that in code defeats the purpose of the function).
arguments
    pv.nrWorkers = str2double(getenv("NS_PARPOOL"));
    pv.spmdEnabled = false
end
pool = gcp("nocreate");
if isnan(pv.nrWorkers)
    % The environment variable was not set: use a pool only if it already 
    % exists     
elseif pv.nrWorkers==0
    % Explicitly disallow using a pool (but do not close it).
    pool = [];
elseif isempty(pool)
    % create a pool with the requested number of workers
    pool = parpool('Processes',pv.nrWorkers,'SpmdEnabled',pv.spmdEnabled);
else 
    % use the pool as is. 
end       

end
