function [pv,ok] = setupPython(env)
% Resolve and, when needed, load Python from the configured environment.
%
% NS_PY_ROOT specifies the folder that directly contains the Conda/Mamba/
% micromamba environments. NS_PY_ENV specifies the default environment name. An environment
% name can also be passed explicitly to this function. The first output is a
% structure describing the environment; the second output reports success.
%
% On Unix python is set to OutOfProcess (the only way to avoid library
% conflicts and segfaults in my experience)
arguments
    env (1,1) string = string(getenv("NS_PY_ENV"))
end

assert(isempty(getenv("NS_CONDA")), ...
    "NS_CONDA is obsolete. Set NS_PY_ROOT to the folder containing the environments and NS_PY_ENV to the environment name.");
assert(isempty(getenv("NS_PYTHON")), ...
    "NS_PYTHON is obsolete. Set NS_PY_ROOT to the folder containing the environments and NS_PY_ENV to the environment name.");

root = string(getenv("NS_PY_ROOT"));
assert(strlength(root) > 0, ...
    "NS_PY_ROOT is not set. Set it to the folder directly containing the Conda/Mamba/micromamba environments.");
assert(exist(root,"dir") == 7, ...
    "NS_PY_ROOT does not exist: %s",root);
assert(strlength(env) > 0, ...
    "NS_PY_ENV is not set. Set it to the environment name, or pass an environment name to setupPython.");

prefix = fullfile(root,env);
assert(exist(prefix,"dir") == 7, ...
    "Python environment %s does not exist under NS_PY_ROOT: %s",env,prefix);
if ispc
    p.executable = string(fullfile(prefix,"python.exe"));
    p.bin = string(fullfile(prefix,"Scripts"));
    p.extension = ".exe";
else
    p.executable = string(fullfile(prefix,"bin","python"));
    p.bin = string(fullfile(prefix,"bin"));
    p.extension = "";
end
assert(exist(p.executable,"file") == 2, ...
    "Python executable not found for environment %s: %s",env,p.executable);
p.root = root;
p.name = env;
p.prefix = string(prefix);

pv = pyenv;
if ~strcmpi(string(pv.Executable),p.executable)
    if pv.Status == "Loaded" && pv.ExecutionMode == "OutOfProcess"
        terminate(pyenv);
        pv = pyenv;
    end

    if pv.Status == "NotLoaded"
        if isunix
            pyenv('ExecutionMode','OutOfProcess');
        end
        pyenv('Version',p.executable);
    else
        error('Could not select Python environment %s.',env);
    end
end

pv = pyenv();
assert(strcmpi(string(pv.Executable),p.executable), ...
    "Python setup selected %s instead of the requested environment executable %s.",pv.Executable,p.executable);
fprintf('Python setup successful (Exe: %s Ver:%s)\n',pv.Executable,pv.Version)
ok = pv.Version ~= "";


