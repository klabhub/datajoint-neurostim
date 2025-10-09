function [signal,time,channelInfo,recordingInfo] = cascade(key,parms)
% Spike inference from calcium data using the CASCADE deep neural network.
%
% Rupprecht, P., Carta, S., Hoffmann, A. et al.
% A database and deep learning toolbox for noise-optimized, generalized
% spike inference from calcium imaging.
% Nat Neurosci 24, 1324–1337 (2021).
% https://doi.org/10.1038/s41593-021-00895-5
%
% SETUP
% % 1. Install micromamba (or miniconda, miniforge, but micromamba seems to work
%   better on HPC)
% 2. Create the cascade environment  (this is the non-gpu version and uses
% the defaults channel because these old TF versions are not available on
% the conda-forge)
% micromamba create -n cascade -c defaults --strict-channel-priority
%   python=3.7 "tensorflow=2.3.*" "keras=2.3.1"   "h5py<3" "numpy<1.20"
%   "scipy<1.6" matplotlib seaborn ruamel.yaml spyder
%
% 3.  Set your NS_PYTHON environment variable to the executable that can run
% python in an enviroment.  For instance '~/miniforge3/bin/conda' or
% '~/.local/bin/micromamba'
% This will be used to construct a bash command that calls python in the
% cascade enviroment.
%
% 4. Clone the Cascade repo to your home folder:
%  git clone  https://github.com/HelmchenLabSoftware/Cascade.git
% Store this folder in the NS_CASCADE environment variable. If this is not
% defined, Matlab will assume that this repo is in the home folder (i.e.
% ~/Cascade)
%
% USAGE
% Define a row in ns.CParm with the following parms struct
%
% parms.model = 'Global_EXC_15Hz_smoothing100ms';
% See https://github.com/HelmchenLabSoftware/Cascade/blob/master/Pretrained_models/available_models.yaml
% for a list of available models.
% parms.baseline = 'maximin' % How to determine F0
% parms.sigma = 10; % In seconds ; Gaussian sigma for F0
% parms.window = 60 ; % in seconds; window for maximim
% parms.count    = false; % Predict spike counts instead of spike  probability
% parms.fluorescence = 'fluorescence' ;  % The C channel with the neuropil corrected, raw  F
% And add this to the CParm table with a call like this
% cascade = struct('ctag','cascade',...       % Name of this C
%                  'description','Cascade inferred spikes',...
%                  'extension','.sbx',...  % Files to process
%                  'fun','sbx.cascade',...      % Use this function to do the work
%                  'parms',parms);
% insertIfNew(ns.CParm,cascade);
%
% This instructs populate(ns.C,'ctag="cascade"') to fill the ns.C table
% with spiking probabilities (count=false) for each sbx file. If such file
% does not have a "fluorescence" row in ns.C, this will fail.
% The typical workflow is to populate sbx.Preprocessed first (e.g. using
% suite2p), then the fluorescence rows in ns.C (neuropil corrected F, see
% sbx.read), and then do spike inference to fill the "cascade" rows of
% ns.C.
% parms.options  is  an optional struct with parameters that will be passed to the
% cascade.py script:
% parms.options.verbosity  (defaults to 0)
% parms.options.model_folder  (defaults to 'Pretrained_models')
% parms.options.threshold (defaults to 0)
%
% BK - Oct 2025
arguments
    key (1,1) struct
    parms (1,:) struct
end

assert(all(ismember(["model" "baseline" "sigma" "window" "fluorescence"],fieldnames(parms))),"The parms struct is missing essential elements.")
dffFile = tempname + ".mat"; % This will be a temp save file for dFF
resultsFile = strrep(dffFile,".mat","_cascade.mat"); % File where the cascade python code saves the results
% Query the CChannel table for fluorescence time series
allF = (ns.C & struct('ctag',parms.fluorescence) & rmfield(key,'ctag') )* (ns.CChannel &  proj(sbx.PreprocessedRoi,'roi->channel'));
nrRoi = count(allF);
assert(nrRoi>0,"No fluorescence data (ctag=%s) found for %s on %s at %s",parms.fluorescence,key.subject,key.session_date,key.starttime)
% Check that cascade is installed
cascadeFolder = getenv("NS_CASCADE");
if isempty(cascadeFolder)
    if ispc
        cascadeFolder =fullfile(getenv('USERPROFILE') , "Cascade");
    else
        cascadeFolder = fullfile(getenv('HOME'),"Cascade");
    end
end
assert(exist(cascadeFolder,"dir"),"Please install the Cascade github repository or set the NS_CASCADE environment variable (not found at  %s)",cascadeFolder);
% Check that we can run PYTHON in a conda/mamba environment
pythonRunner = getenv("NS_PYTHON");
if isempty(pythonRunner)
    if ispc
        pythonRunner = fullfile(getenv('USERPROFILE') , "miniconda3/condabin/conda.bat");
    else
        pythonRunner = fulllfile(getenv('HOME'),"/miniconda3/condabin/conda");
    end
end
assert(exist(pythonRunner,"file"),"Set the NS_PYTHON environment variable to a conda/mamba/micromamba executable (not found at  %s)",pythonRunner);


%% Get the neuropil corrected fluorescence
% TODO : memory check and loop if not enough memory available.
tpl = fetch(allF);
F = fetch(allF,'signal');
F =cat(2,F.signal);
prep = fetch(sbx.Preprocessed & key,'framerate','nrplanes');
time =  fetch1(allF,'time','LIMIT 1'); %[start stop nrSamples]

%% Determine dFF
% Replacing F to save a bit of memory
switch upper(parms.baseline)
    case "MAXIMIN"
        [~,F]= sbx.baseline_maximin(F,prep.framerate/prep.nrplanes,parms.sigma,parms.window);
    otherwise
        error('Not implemented yet')
end

% Calculate the noise level using the Cascade formula:
% the median absolute dF/F difference between two subsequent time points.
% This is a outlier-robust measurement that converges to the simple
% standard deviation of the dF/F trace for  uncorrelated and outlier-free dF/F traces.
% The value is divided by the square root of the frame rate to make it comparable across recordings with different frame rates.

noiseLevel = 100*median(abs(diff(F,1,1)),1,"omitmissing")./sqrt(prep.framerate);


%% Construct the python command.
cfd = fileparts(mfilename('fullpath'));
toolsPath = strrep(fullfile(fileparts(cfd),'tools'),"\","/");
cascadeFolder = strrep(cascadeFolder,"\","/");
pyCmd = sprintf('"%s/cascade.py" "%s" "%s" --cascade_folder "%s"',toolsPath,dffFile,parms.model,cascadeFolder);
if parms.count 
    % Only add this step if we're computing the count (and not the
    % probability)
    pyCmd = pyCmd + " --discrete_spikes ";
end
% add options
if isfield(parms,'options')
    fn = fieldnames(parms.options);
    opts  ="";
    for i=1:numel(fn)
        opts = opts + " --" + fn{i} + " " + string(parms.options.(fn{i}));
    end
    pyCmd =pyCmd + opts;
end
pyEnv = pythonRunner + " run -n cascade python ";
cmd = sprintf('%s %s',pyEnv,pyCmd);

%% Run cascade in chunks
nrSamples  = time(3);
samples= 0:nrSamples;
signal = zeros(nrSamples,nrRoi);
CHUNK = 200; % With 70k samples, this results in <20GB of Ram and ~ 20 cores
for ix =1:CHUNK:nrRoi
    thisRoi = ix:min(ix+CHUNK-1,nrRoi);
    thisNrRoi = numel(thisRoi);
    dFF = F(:,thisRoi);
    save(dffFile,"dFF");
    % Execute
    fprintf('Starting Cascade with %s on %d ROIs with %d samples @%s\n',parms.model,thisNrRoi,nrSamples,datetime('now'));
    [status,msg] = system(cmd ,'-echo');
    assert(status==0,'Cascade failed with message %s',msg);
    fprintf('Cascade completed successfully.\n');
    %% Read the results file
    % The cascade.py script in the tools folder puts its output in a file with
    % the same name as the dffFile but _cascade.mat suffix.
    % It containst pSpike ; a matrix with spiking probabilities for each time
    % point (row) and neuron (col) - matching F. And sSpike; a cell array [1
    % nrNeurons] with each cell the samples in F at which a spike occurred.
    assert(exist(resultsFile,"file"),"The Cascade output file %s could not be found.",resultsFile);
    if parms.count
        % Convert spike times into a time course of counts.
        sSpike  = load(resultsFile,'sSpike').sSpike;
        assert(thisNrRoi == numel(sSpike),'Missing roi from cascade results?')
        for j=1:numel(sSpike)
            signal(:,thisRoi(j)) = groupcounts(sSpike{j},samples,IncludeEmptyGroups=true);
        end
        name = 'spikeCount';
    else
        % Return the spiking probability as the signal
        load(resultsFile,'pSpike');
        assert(thisNrRoi == size(pSpike,2),'Missing roi from cascade results?')
        signal(:,thisRoi) = pSpike;
        name = 'spikeProbability';
    end
end
channelInfo =struct('name',name,'nr',{tpl.channel}','noiseLevel',num2cell(noiseLevel)');
recordingInfo =struct('source','cascade');
%% Clean up temp files.
delete(dffFile);
delete(resultsFile);

