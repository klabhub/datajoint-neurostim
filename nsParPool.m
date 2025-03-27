function pool = nsParPool()
% Function to start a local parpool based on an environment variable. 
% 
% This is useful to avoid storing things like nrWorkers with parameters for
% a makeTuples function  (which forces the same number of workers forever
% and is annoying when debugging).
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
% If you want to use the parfor with 10 workers, call setenv("NS_PARFOR",10) 
% or you start a parpool directly with a call to parpool() and ensure that
% NS_PARFOR environment variable is not set (or "")
%
% To temporarily use the for-loop even with a parpool running, call
% setenv("NS_PARFOR",0)
%
nrWorkers = str2double(getenv("NS_PARFOR"));
pool = gcp("nocreate");
if isnan(nrWorkers)
    % The environment variable was not set: use a pool only if it already 
    % exists     
elseif nrWorkers==0
    % Explicitly disallow using a pool (but do not close it).
    pool = [];
elseif isempty(pool)
    % create a pool with the requested number of workers
    pool = parpool('Processes',nrWorkers);
else 
    % use the pool as is. 
end       

end