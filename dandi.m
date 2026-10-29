function dandi(dandiSet, root, pv)
% Validate and upload a folder to the Dandi archive
% Run this after ns.Experiment/nwbExport. 
% 
% Dandi documentation https://github.com/dandi/dandi-docs
% 
% EXAMPLE
% dandi (00170,"c:/temp/dandi",pythonEnvironment ="dandi",staging=false)
%  This will validate, organize and upload the data in c:/temp/dandi/00170
%  to the Dandi archive.
% 
% This function wil first call dandi download, then organize, then upload.
% Each of those calls can be modified by specifying downloadOptions,
% uploadOptions, or organizeOptions as strings. See dandi --help for details, but for
% instance:
% The uploadOptions can be set to 
% "--sync" to  Delete assets on the server that do not
%                                  exist locally 
% --existing OVERWRITE : force upload
%                                  if either size or modification time differs
% --existing REFRESH : upload only if local
%                                  modification time is ahead of the remote.
%
% This uses the environment configured by NS_PY_ROOT and NS_PY_ENV (or the
% explicit pythonEnvironment parameter) and requires NWB/Dandi to be installed.
% installed. 
% I used this to install it:
% conda create -n nwb python=3.10
% conda activate nwb 
% pip install nwbinspector
% pip install "nwbinspector==0.6.3" "pynwb>=3.0" "dandi>=0.60"
%
% This worked fine on a Windows/miniconda install, but on an HPC system it
% failed due to a missing deno installation. There I had to install
% dandi=0.68.1 to get it to work.
%
% Note that just installing default dandi (using conda install) led to
% run time errors due to a missing hdmf.array class. Not sure why.
%
% SEE ALSO ns.Experiment/nwbExport
arguments
    dandiSet (1,1) string   % The Dandi Set unique identifier as a string
    root (1,1) string       % The folder that contains a folder called export with the nwb data files exported by ns.Experiment/nwbExport
    pv.pythonEnvironment (1,1) string = string(getenv("NS_PY_ENV"))
    pv.staging (1,1) logical = false    % Set to true if this is a staging DandiSet
    pv.upload (1,1) logical = true     % Set to false to skip the upload step 
    pv.uploadOptions (1,1) string = ""   % Provide additional verbatim input to the dandi upload commands
    pv.downloadOptions (1,1) string = "--existing refresh" % Provide additional verbatim input to the dandi download commands
    pv.organizeOptions (1,1) string = "--files-mode auto" % Provide additional verbatim input to the dandi organize commands
    pv.runNwbInspector (1,1) logical = true; % Set to false to skip the initial validation step. (Validation is always run after dandi organize)
end
if pv.staging 
    url = 'api-staging.dandiarchive.org/api';
else
    url = 'dandiarchive.org';
end


% The data are in root/export
exportFolder = fullfile(root,"export");
assert(exist(exportFolder,"dir"),"%s does not exist. Run nwbExport(ns.Experiment)?",exportFolder)
%% Run nwbinspector dandi validation in a conda environment
if pv.runNwbInspector
fprintf('**** Running nwbinspector...\n');                
    runInPython("nwbinspector --config dandi .", exportFolder,pv.pythonEnvironment);
end
%% Delete the existing dandiset folder (otherwise organize generates errors)
dandiSetFolder = fullfile(root,dandiSet);
if exist(dandiSetFolder,"dir")
    rmdir(dandiSetFolder,'s');
end
%% Download .yaml /data to prepare for upload. 
% This wil create the dandi set folder root\dandiSetNumber with 
% the dandiset.yaml file in it.
fprintf('**** Downloading dandiset ...\n');                    
dandiCmd = sprintf("dandi download " + pv.downloadOptions + " --output-dir .. --download dandiset.yaml https://%s/dandiset/%s/draft/",url,dandiSet);
runInPython(dandiCmd,exportFolder,pv.pythonEnvironment);
    
%% Organize the files and revalidate
% Force session_id in the name of the file as otherwise
% dandi tries to rename the file to subject.nwb and
% when later sessions are added, they get random
% suffixes to disambiguate. This seems more readable.
% Note that session_id includes the starttime (i.e. the file corresponds to
% an experiment in the Neurostim sense).
% This will create symlinks or copies renamed according to the BIDS/Dandi
% conventions
fprintf('**** Organizing and validating dandiset %s ...\n',dandiSet);                

dandiCmd = sprintf("dandi organize " + pv.organizeOptions + " --dandiset-path ../%s --required-field session_id %s",dandiSet,exportFolder);
runInPython(dandiCmd,exportFolder,pv.pythonEnvironment);    

         
%% Upload
fprintf('**** Uploading dandiset %s ...\n',dandiSet); tic;
% On dandi 0.60 I had to add these command line arguments 

if pv.staging
    dandiCmd = "dandi upload " + pv.uploadOptions  + " -i dandi-staging .";          
else
    dandiCmd = "dandi upload " + pv.uploadOptions  +  " -i dandi .";   
end

runInPython(dandiCmd,dandiSetFolder,pv.pythonEnvironment);
fprintf('*** Upload complete (%s)\n',seconds(toc));
end

%% Helper function to run a dandi command in a specified folder using a 
% user defined conda environment.
function status = runInPython(dandiCmd,wrkFldr,pythonEnvironment)
arguments
    dandiCmd (1,1) string  % Command to run
    wrkFldr (1,1)  string   % Working folder for command    
    pythonEnvironment (1,1) string = string(getenv("NS_PY_ENV")) % Python environment to use
end 
py = setupPython(pythonEnvironment);
tokens = regexp(char(dandiCmd),'^(?<name>\S+)(?<args>.*)$','names');
tool = fullfile(py.bin,string(tokens.name) + py.extension);
assert(exist(tool,"file") == 2, ...
    "Python command %s not found in environment %s: %s",tokens.name,py.name,tool);
if ispc
    cmd = sprintf('cd /d "%s" & "%s"%s',wrkFldr,tool,tokens.args);
else
    cmd = sprintf('cd "%s" && "%s"%s',wrkFldr,tool,tokens.args);
end

status =  system(cmd,'-echo');
assert(status==0,"%s failed",dandiCmd);

end
