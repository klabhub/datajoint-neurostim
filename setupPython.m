function ok = setupPython
% Check that Python is accessible and try to connect if not.
% Order:
% 1. The python that is already loaded (pyenv)
% 2. The python in the conda environment 'matlab'
% 3. The python in the default conda environment 
%  
% The conda install is determined from the NS_CONDA environment variable.
% 
pv = pyenv;
if ispc
    exe = 'python.exe';
else
    exe = 'python';
end
if isempty(pv.Version)
    conda = getenv("NS_CONDA");
    if ~isempty(conda) && exist(conda,'dir')
        % Check for a conda environment called matlab.
        f = fullfile(conda,'envs','matlab','bin',exe);
        if exist(f,'file')
            pv = pyenv('Version',f);
        else
            % Use the default install of miniconda
            f= fullfile(conda,'bin',exe);
            pv = pyenv('Version',f);
        end
    else
        fprintf('CONDA dir %s does not exist. Define NS_CONDA.',conda);
        f ='';
    end
end

if isempty(pv.Version)
    fprintf(2,'Python setup unsuccessful (conda: %s file: %s)\n',conda,f)
else 
    fprintf('Python setup successful (Exe: %s Ver:%s)\n',pv.Executable,pv.Version)
end
ok = pv.Version ~="";