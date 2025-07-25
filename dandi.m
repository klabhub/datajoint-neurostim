function dandi(dandiSet, root, pv)
% Validate and upload a folder to the Dandi archive
% Run this manually after nwbExport, or let nwbExport call it automatically
% once it has created the .nwb files. 
% 
% Dandi documentation https://github.com/dandi/dandi-docs
% 
% EXAMPLE
% dandi (00170,"c:/temp/dandi",condaEnvironment ="dandi",staging=false)
%  This will validate, organize and upload the data in c:/temp/dandi/00170
%  to the Dandi archive.
% 
% The uploadOptions can be set to 
% --sync to  Delete assets on the server that do not
%                                  exist locally 
% --existing OVERWRITE : force upload
%                                  if either size or modification time differs
% --existing REFRESH : upload only if local
%                                  modification time is ahead of the remote.
%
% This requires a CONDA installation defined by the NS_CONDA environment
% variable (e.g. /home/joe/miniconda3) and a CONDA environment with NWB/Dandi 
% installed. 
% I used this to install it:
% conda create -n nwb python=3.10
% conda activate nwb 
% pip install nwbinspector
% pip install "nwbinspector==0.6.3" "pynwb>=3.0" "dandi>=0.60"
%
% This worked fine on a Windows/miniconda install, but on an HPC system it
% failed due to a missing deno installation. There I had to install
% dandi=0.60 to get it to work (this means no zarr files or other more
% recent improvements to dandi/nwb.)
% Note that just installing default dandi (using conda install) led to
% run time errors due to a missing hdmf.array class. Not sure why.
%
% SEE ALSO ns.Experiment/nwbExport
arguments
    dandiSet (1,1) string   % The Dandi Set unique identifier as a string
    root (1,1) string       % The folder that contains a folder called export with the nwb data files exported by ns.Experiment/nwbExport
    pv.condaEnvironment (1,1) string = "nwb"
    pv.staging (1,1) logical = false    % Set to true if this is a staging DandiSet
    pv.upload (1,1) logical = true     % Set to false to skip the upload step 
    pv.uploadOptions (1,1) string = ""   % Provide additional verbatim input to the dandi upload commands
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
runInConda("nwbinspector --config dandi .", exportFolder,pv.condaEnvironment);
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
dandiCmd = sprintf("dandi download --existing refresh --output-dir .. --download dandiset.yaml https://%s/dandiset/%s/draft/",url,dandiSet);
runInConda(dandiCmd,exportFolder,pv.condaEnvironment);
    
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

dandiCmd = sprintf("dandi organize --files-mode auto --dandiset-path ../%s --required-field session_id %s",dandiSet,exportFolder);
runInConda(dandiCmd,exportFolder,pv.condaEnvironment);    

         
%% Upload
fprintf('**** Uploading dandiset %s ...\n',dandiSet); tic;
% On dandi 0.60 I had to add these command line arguments 
v60Options = "--format PYOUT --path-type EXACT";
if pv.staging
    dandiCmd = "dandi upload " +v60Options + pv.uploadOptions  + " -i dandi-staging .";          
else
    dandiCmd = "dandi upload " +v60Options + pv.uploadOptions  +  " -i dandi .";   
end

runInConda(dandiCmd,dandiSetFolder,pv.condaEnvironment);
fprintf('*** Upload complete (%s)\n',seconds(toc));
end

%% Helper function to run a dandi command in a specified folder using a 
% user defined conda environment.
function status = runInConda(dandiCmd,wrkFldr,condaEnvironment)
arguments
    dandiCmd (1,1) string  % Command to run
    wrkFldr (1,1)  string   % Working folder for command    
    condaEnvironment (1,1) string = "nwb" % Conda environment to use    
end 
condaFldr = getenv('NS_CONDA');
assert(exist(condaFldr,"dir"),'Set the NS_CONDA variable (%s) to point to your Conda installation (e.g. /home/user/miniconda3',condaFldr)
               
if ispc
    cmd = sprintf('"%s\\Scripts\\activate.bat" "%s" & conda activate %s & cd "%s" & %s',condaFldr,condaFldr,condaEnvironment,wrkFldr,dandiCmd);                    
else
    cmd = sprintf(['__conda_setup="$(''%s/bin/conda'' shell.bash hook 2>/dev/null)"; ' ...
    'if [ $? -eq 0 ]; then eval "$__conda_setup"; ' ...
    'elif [ -f "%s/etc/profile.d/conda.sh" ]; then source "%s/etc/profile.d/conda.sh"; ' ...
    'else export PATH="%s/bin:$PATH"; fi; unset __conda_setup; ' ...
    'conda activate %s; cd %s; eval "%s"'], ...
    condaFldr, condaFldr, condaFldr, condaFldr, condaEnvironment, wrkFldr,dandiCmd);
end

status =  system(cmd,'-echo');
assert(status==0,"%s failed",dandiCmd);

end