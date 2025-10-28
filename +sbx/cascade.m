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
% Use a ? wildcard (e.g., parms.model = 'Global_EXC_?Hz_smoothing100ms') to
% use  a model with the frequency that is closest to the actual sampling frequency of the
% data. For instance, this allows automatic adjustment to the 7.5 Hz model
% for 2 plane recordings at 15Hz.
%
% parms.prep = 'suite2p'; % Which sbx.Preprocessed to use
% parms.restrict = 'pcell>0.9'; % Optional restriction on which rois to do
% parms.neuropilFactor = 0.7; % Subtract 0.7*Fneu from F.
% parms.baseline = 'maximin'; % How to determine F0
% parms.sigma = 10; % In seconds ; Gaussian sigma for F0
% parms.window = 60 ; % in seconds; window for maximim
% parms.count    = false; % Predict spike counts instead of spike  probability
% And add this to the CParm table with a call like this
% cascade = struct('ctag','cascade',...       % Name of this C
%                  'description','Cascade inferred spikes',...
%                  'extension','.sbx',...  % Files to process
%                  'fun','sbx.cascade',...      % Use this function to do the work
%                  'parms',parms);
% insertIfNew(ns.CParm,cascade);
%
% This instructs populate(ns.C,'ctag="cascade"') to fill the ns.C table
% with spiking probabilities (count=false) for each sbx file.
% The typical workflow is to populate sbx.Preprocessed first (e.g. using
% suite2p), and then do spike inference to fill the "cascade" rows of
% ns.C.
%
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

%% Prepare to run
prep = fetch(sbx.Preprocessed & key & struct('prep', parms.prep),'*');
assert(~isempty(prep),'No preprocessed data for %s in session %s for subject %s. Run populate(sbx.Preprocessed,prep="%s") first',parms.prep,key.session_date,key.subject,parms.prep);
warning('off','backtrace');
assert(all(isfield(parms,["model" "baseline" "sigma" "window" "count"])),'cascade parameters incomplete.');
if ~isfield(parms,'neuropilFactor')
    parms.neuropilFactor = 0.7; % Default neuropil correction factor
end
if ~isfield(parms,'perExperiment')
    parms.perExperiment = false;
end
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

rate  = (prep.framerate/prep.nrplanes); % Match dt to framerate
if contains(parms.model,"?Hz")
    % Find the model with the closest matching frequency (mainly intended to adjust
    % for two plane recordings that have an effective sampling rate that
    % is half the 15Hz of the microscope)
    models = getCascadeModels(cascadeFolder);
    match = regexp(models,strrep(parms.model,'?Hz','(?<rate>\d+)Hz'),'names');
    availableRate= cellfun(@(x) str2double(x.rate),match(~cellfun(@isempty,match)));
    assert(~isempty(availableRate),'No matching Cascade models for %s',parms.model);
    [~,ix] =min(abs(rate-availableRate));
    parms.model = strrep(parms.model,'?',num2str(availableRate(ix)));
end

match = regexp(parms.model,'_(?<rate>\d+)Hz','names');
assert(~isempty(match),'Cannot extract sampling rate from model name %s',parms.model)
modelRate =str2double(match.rate);
ratio = modelRate/rate;
if ratio >2 || ratio < .5
    warning('Large mismatch between sampling rate %.2f and Cascade model %.2f',rate,modelRate);
end




%% Construct the python command.
cfd = fileparts(mfilename('fullpath'));
toolsPath = strrep(fullfile(fileparts(cfd),'tools'),"\","/");
cascadeFolder = strrep(cascadeFolder,"\","/");
tmpDffFile = tempname + ".mat"; % This will be a temp save file for dFF
pyCmd = sprintf('"%s/cascade.py" "%s" "%s" --cascade_folder "%s"',toolsPath,tmpDffFile,parms.model,cascadeFolder);
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



%% Get the frames for this experiment
[keepFrameIx,frameNsTime] = sbx.framesForExperiment(key);
nrFrames = numel(keepFrameIx);
fldr= fullfile(folder(ns.Experiment & key),fetch1(sbx.Preprocessed & key & struct('prep',parms.prep),'folder'));
allPlanesSignal = [];
allPlanesChannelInfo = [];

for pl=0:prep.nrplanes-1
    planeFolder = fullfile(fldr,"plane" +  string(pl));
    if parms.perExperiment
        cascadeResultsFilename = fullfile(planeFolder ,key.starttime + "." + key.ctag + ".cascade.mat");
    else
        cascadeResultsFilename = fullfile(planeFolder ,key.ctag + ".cascade.mat");
    end
    if exist(cascadeResultsFilename,"file")
        fprintf('Cascade results already exist for plane %d in %s. Reading from file.\n',pl,fldr);
        load(cascadeResultsFilename,'signal','channelInfo');
    else
        if isfield(parms,'restrict')
            restrict = parms.restrict;
        else
            restrict = 'true';
        end
        % Restrict with a query on sbx.PreprocessedRoi
        roi = [fetch((sbx.PreprocessedRoi & restrict & struct('plane',pl)) & key ,'roi').roi]';

        if  isempty(roi)
            fprintf("No ROI match the restriction. Skipping.");
            continue; % Skip to next plane
        else
            roiOffset = fetch1(aggr(sbx.Preprocessed&key,sbx.PreprocessedRoi & key & sprintf('plane<%d',pl),'max(roi)->offset'),'offset');
            roiOffset(isnan(roiOffset)) =0; % Plane 0 -> no offset
            try
                [signal,channelInfo]= runCascade(roi,roiOffset,parms,prep,planeFolder,cascadeResultsFilename,cmd,tmpDffFile,keepFrameIx);
            catch me
                if contains(me.message,"Python process terminated unexpectedly")
                    fprintf('Python terminated unexpectedly. Retrying...')
                    terminate(pyenv);
                    [signal,channelInfo]= runCascade(roi,roiOffset,parms,prep,planeFolder,cascadeResultsFilename,cmd,tmpDffFile,keepFrameIx);
                end
            end
            fprintf('Completed cascade inference on %d rois. \n',size(signal,2));
        end
    end

    if parms.perExperiment
        % results file only contains the relevant frames - no selection
        % needed.
    else
        % Select the frames for this experiment from the signal
        framesNotInFile = sum(keepFrameIx > size(signal,1)); %#ok<*UNRCH>
        assert(framesNotInFile==0,"%s has %d too few frames for %s on %s",cascadeResultsFilename,framesNotInFile,key.starttime, key.session_date);
        signal = signal(keepFrameIx,:);
    end
    % Concatenate across planes
    allPlanesSignal = [allPlanesSignal  signal]; %#ok<AGROW>
    allPlanesChannelInfo = [allPlanesChannelInfo  ;channelInfo]; %#ok<AGROW>
end

% Rename for output values
signal = allPlanesSignal;
channelInfo = allPlanesChannelInfo;
time = [frameNsTime(1) frameNsTime(end) nrFrames];
recordingInfo = struct('dummy',true);

end

function [signal,channelInfo] = runCascade(roi,roiOffset,parms,prep,planeFolder,cascadeResultsFilename,cmd,tmpDffFile,keepFrameIx)
% Run cascade on the fluorescence data in planeFolder for the rois in keepRoi
% and save the results in cascadeResultsFilename
arguments
    roi(:,1) double % roi numbers (across planes) to process
    roiOffset (1,1) double % Offset for the current plane
    parms (1,1) struct
    prep (1,1) struct
    planeFolder (1,1) string
    cascadeResultsFilename (1,1) string
    cmd (1,1) string
    tmpDffFile (1,1) string
    keepFrameIx (1,:) double; % Only used when parms.perExperiment =true
end

%% Load the F data and subtract neuropil
thisFile = fullfile(planeFolder,'F.npy');
if ~exist(thisFile,"file")
    error('File %s does not exist',thisFile);
end
F =  ndarrayToArray(py.numpy.load(thisFile,allow_pickle=true),single=true);
roiIx = roi-roiOffset;
F = F(roiIx,:)';
thisFile = fullfile(planeFolder,'Fneu.npy');
if ~exist(thisFile,"file")
    error('File %s does not exist',thisFile);
end
Fneu =  ndarrayToArray(py.numpy.load(thisFile,allow_pickle=true),single=true);
Fneu = Fneu(roiIx,:)';
F = F-parms.neuropilFactor*Fneu;

if parms.perExperiment
    % Per Experiment mode
    % Select the frames for this experiment from the signal
    framesNotInFile = sum(keepFrameIx > size(F,1)); %#ok<*UNRCH>
    assert(framesNotInFile==0,"%s has %d too few frames for %s on %s",cascadeResultsFilename,framesNotInFile);
    F= F(keepFrameIx,:);
end

%% Determine dFF
% Replacing F to save a bit of memory
switch upper(parms.baseline)
    case "MAXIMIN"
        [F0,F]= sbx.baseline_maximin(F,prep.framerate/prep.nrplanes,parms.sigma,parms.window);
    otherwise
        error('Not implemented yet')
end

% Calculate the noise level using the Cascade formula:
% the median absolute dF/F difference between two subsequent time points.
% This is a outlier-robust measurement that converges to the simple
% standard deviation of the dF/F trace for  uncorrelated and outlier-free dF/F traces.
% The value is divided by the square root of the frame rate to make it comparable across recordings with different frame rates.

noiseLevel = 100*median(abs(diff(F,1,1)),1,"omitmissing")./sqrt(prep.framerate/prep.nrplanes);
%  Use the roi number in the channelInfo to keep consistent naming with
%  PreprocessedRoi
channelInfo =  struct('nr',num2cell(roi),'noiseLevel',num2cell(noiseLevel'),'model',char(parms.model));

%% Run cascade in chunks
nrRoi = size(F,2);
samples= 0:prep.nrframesinsession;
signal = [];
CHUNK = 200; % With 70k samples, this results in <20GB of Ram and ~ 20 cores
tmpResultsFile = strrep(tmpDffFile,".mat","_cascade.mat"); % File where the cascade python code saves the results
for ix =1:CHUNK:nrRoi
    thisRoi = ix:min(ix+CHUNK-1,nrRoi);
    thisNrRoi = numel(thisRoi);
    dFF = F(:,thisRoi);
    save(tmpDffFile,"dFF");
    % Execute
    fprintf('Starting Cascade with %s on %d ROIs with %d samples @%s\n',parms.model,thisNrRoi,prep.nrframesinsession,datetime('now'));
    [status,msg] = system(cmd ,'-echo');
    assert(status==0,'Cascade failed with message %s',msg);
    fprintf('Cascade completed successfully.\n');
    %% Read the results file
    % The cascade.py script in the tools folder puts its output in a file with
    % the same name as the dffFile but _cascade.mat suffix.
    % It containst pSpike ; a matrix with spiking probabilities for each time
    % point (row) and neuron (col) - matching F. And sSpike; a cell array [1
    % nrNeurons] with each cell the samples in F at which a spike occurred.
    assert(exist(tmpResultsFile,"file"),"The Cascade output file %s could not be found.",tmpResultsFile);
    if parms.count
        % Convert spike times into a time course of counts.
        sSpike  = load(tmpResultsFile,'sSpike').sSpike;
        assert(thisNrRoi == numel(sSpike),'Missing roi from cascade results?')
        for j=1:numel(sSpike)
            signal  = [signal  groupcounts(sSpike{j},samples,IncludeEmptyGroups=true)]; %#ok<AGROW>
        end
    else
        % Use the spiking probability as the signal
        load(tmpResultsFile,'pSpike');
        assert(thisNrRoi == size(pSpike,2),'Missing roi from cascade results?')
        signal = [signal pSpike]; %#ok<AGROW>
    end
end

%% Save results for future reference
save(cascadeResultsFilename, 'signal','channelInfo','F0');
fprintf('Saved results to %s\n',cascadeResultsFilename);

%% Clean up temp files.
delete(tmpDffFile);
delete(tmpResultsFile);
end



function modelNames = getCascadeModels(cascadePath)
% Read available CASCADE models from YAML file

filePath = fullfile(cascadePath,'Pretrained_models','available_models.yaml');

% Read file
fileID = fopen(filePath, 'r');
fileText = fread(fileID, '*char')';
fclose(fileID);

% Extract model names
lines = splitlines(fileText);
modelNames = {};

for i = 1:length(lines)
    line = strtrim(lines{i});
    % Model names are at the start of line and end with ':'
    if ~isempty(line) && ~startsWith(line, '#') && ...
            ~startsWith(line, ' ') && endsWith(line, ':')
        modelName = strrep(line, ':', '');
        modelNames{end+1} = modelName; %#ok<AGROW>
    end
end
end


