function ok = setupPython(env)
% Check that Python is accessible and try to connect if not.
% Order:
% 1. The python that is already loaded (pyenv)
% 2. The python in the conda environment 'matlab'
% 3. The python in the conda environment specified as the input argument
% 4. The base conda environemnt (if 2 or 3 do not exist)
%
% The conda install is determined from the NS_CONDA environment variable.
% Set it to something like "~/miniconda3"
arguments
    env (1,1) string = "matlab" % Default env, to be created by user
end
pv = pyenv;
if ispc
    exe = 'python.exe';
else
    exe = 'bin/python';
end

envLoaded = extractAfter(pv.Home,'envs/');
if envLoaded ~=env
    % Change of environemtn
    if pv.Status == "Loaded"  && pv.ExecutionMode =="OutOfProcess"
        terminate(pyenv);  % Terminate the old
        pv= pyenv;
    end

    if pv.Status == "NotLoaded"
        % Either just terminated (outtofprocess) or never loaded
        % (inprocess) - set up new
        conda = getenv("NS_CONDA");
        if ~isempty(conda) && exist(conda,'dir')
            % Check for a conda environment called matlab.
            f = fullfile(conda,'envs',env,'bin',exe);
            if exist(f,'file')
                pyenv('Version',f);
            else
                fprintf('CONDA env %s does not exist.',env);
                % Use the default install of miniconda
                f= fullfile(conda,'bin',exe);
                pyenv('Version',f);
            end
        else
            fprintf('CONDA dir %s does not exist. Define NS_CONDA.',conda);
            f ='';
        end
    end
end

% Check the outcome
pv =pyenv();
envLoaded = extractAfter(pv.Home,'envs/');
if envLoaded ~= env
    fprintf(2,'%s could not be loaded, instead we have %s. Hoping for the best.\n',env,envLoaded);
end

if isunix
    % Crashing without this
    pyenv('ExecutionMode','OutOfProcess');
end

if isempty(pv.Version)
    fprintf(2,'Python setup unsuccessful (conda: %s file: %s, env , %s)\n',conda,f,env)
else
    fprintf('Python setup successful (Exe: %s Ver:%s)\n',pv.Executable,pv.Version)
end
ok = pv.Version ~="";