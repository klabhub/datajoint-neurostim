function runsuite2p(opts,db,pv)
% Runs Suite2p using Python. 
% The recommended usage is to specify a condaEnv in which deeplabcut is
% installed, optionally with some bash code to initialize conda. 
%
% Alternatives are to use a singularity container or to just execute the
% python code on the default command line.
arguments
    opts 
    db 
    pv.singularity (1,1) string =""         % THe name of the singularity container to use
    pv.condaEnv (1,1) string =""            % The name of the conda environment to use
    pv.condaInit (1,1) string =""            % A (bash) command to execute before activating the conda environment
end
% command a set of  parms defining how to run

if canUseGPU
    % In case we really need to determine which gpu is availabel, something like this may work
    %  [~,ps] = system('ps -u bart');
    %   ps =strsplit(ps ,{' ','\n'});
    %   ps(cellfun(@isempty,ps)) =[];
    %   ps = reshape((ps),4,[])';
    %   T= table(ps(2:end,1),ps(2:end,2),ps(2:end,3),ps(2:end,4),'VariableNames',ps(1,:));
    % And compare that with
    % nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory --format=csv    
    nvoption = '--nv'; % Singularity only
else    
    nvoption = '';
end

% Save the dicts to tempfiles
optsFile= tempname;
dbFile= tempname;
py.run('import json')
py.run(' with open(filename, ''w'') as file: json.dump(dictionary, file)',dictionary =opts,filename= optsFile)
py.run(' with open(filename, ''w'') as file: json.dump(dictionary, file)',dictionary =db,filename= dbFile)

% Construct the Python command
pythonCmd = sprintf(['from  suite2p import run_s2p',...
          'print(''Reading ops and db, then starting suite2p'');', ...
          'ops=load_dict_from_file(%s);',...
           'db =load_dict_from_file(%s);',...
           'print(ops);',...
            'print(db);',...
            'run_s2p(ops, db);'],optsFile,dbFile);

if pv.singularity ~=""
    % Runin a singularity container
    cmd = sprintf('singularity exec %s %s python -Wdefault -c "%s"',nvoption,parms.singularity, pythonCmd);
elseif pv.condaEnv ~=""
    % Run in a conda environment 
    cmd = sprintf('conda activate %s; python -Wdefault -c "%s"', pv.condaEnv,pythonCmd);
    if pv.condaInit ~=""
        % Prepend conda initialization code provided in the
        % parms
        cmd = pv.condaInit + ";" + cmd;
    end
else
    % Python install without an environment.
    cmd = sprintf('python -Wdefault -c "%s"', pythonCmd);
end

%% Start 
try
    fprintf('Runing suite2p :\n\n %s \n\n',cmd);
    [status] = system(cmd,'-echo');
catch me
    status = -1;
    fprintf('Suite2p failed: %s\n',me.message);
end

if status~=0
    fprintf('Suite2p returede status %d\n',status);
end
end