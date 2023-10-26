function[signal,time,channelInfo,recordingInfo] = read(key,parms)
% Read routine for ns.Continuous to read from Eyelink edf data files, 
% and preprocess them according to the parameters passed in the parms struct, 
% which should have the following fields:
% 
% .channel = cell array with fieldnames that map to eye position/pupil size
% parameters. For instance {'px','py','pa'} will read the raw pupil
% measures from the EDF file and return them as signal columns 1,2,3.
% .eye  = Which eye to use 1=left,2=right
% .downsample = frequency target for downsampling using decimate (optional)
% .fillmissing = Arguments to replace NaN. Passed to fillmissing (e.g. use
% {"nearest","EndValues","nearest"} to replace all NaN with the nearest
% sample.
% 
% See Also fillmissing, decimate, ns.Continuous
arguments
    key % The key of the Expriment table (File and ContinuousParm tuple)
    parms (1,1) struct  =struct% The preprocessing parameters
end

%% Determine which file to read for this experiment
import ns.eyelink.*

qry = ns.File & key;
nrFiles = count(qry);
if nrFiles ~=1
    % Zero or more than 1 file
    error('This experiment has %d files. Cannot proceed.',nrFiles);
else
    % Fetch the file to read
    filename = fullfile(folder(ns.Experiment &key),fetch1(qry,'filename'));
end
if exist(filename,'file') && ~exist(filename,'dir')
    fprintf('Reading from  %s\n',filename);
else
    error('Eyelink %s does not exist',filename);
end


%% The file exist. Open it and read using a MEX file from Christopher Kovach
% Some of the Eyelink messages evaluate to  a matlab function (ellipse)
% but they are upper case. They are evaluated because that is how matlab
% determines whether it can convert something to a number.  I switch off
% the warning here, and back on after loading.
stts = warning;
warning('off','MATLAB:dispatcher:InexactCaseMatch');
data = ns.eyelink.edfmex(char(filename));
% Error checking
notRecorded = ~isfield(data.FSAMPLE,parms.channel);
if any(notRecorded)
   error('Channels %s were not recorded in edf file %s',strjoin(parms.channel(notRecorded),'/'),filename);
end
notRecorded = ~ismember(parms.eye,[data.RECORDINGS.eye]);
if any(notRecorded)
    error('Eye %d were not recorded in edf file %s',parms.eye(notRecorded),filename);
end
warning(stts);

%% Determine alignment between eyelink and neurostim
% Extract TRIALID message
% At the start of each trial; 
% the neurostim eyelink plugin 
% sends the TRIALID message 
% Requests the eyelink clock time and stores this as eyeClockTime
% Hence the eyeClockTime allows us to synchronize the two clocks. 
nrTrials =fetch1(ns.Experiment & key,'trials');
messages= {data.FEVENT.message};
stay = cellfun(@ischar,messages);
messages = messages(stay);
messageEvents = [data.FEVENT(stay)];
match = regexp(messages,'^TRIALID\s+(?<condition>\d+)-(?<trial>\d+)','names');
isTrialID=~cellfun(@isempty,match);
match = [match{:}];
trials= str2num(char(match.trial))'; %#ok<ST2NM>
missingTrials = ~ismember(1:nrTrials,trials);
trialIDEvents = messageEvents(isTrialID);
trialIDTimeEyelink  = nan(nrTrials,1);
trialIDTimeEyelink(trials) = [trialIDEvents.sttime]; % The time of all TRIALID events on the eyelink clock
% Get the equivalent time from Neurostim.
trialIDTimeNeurostim = get(ns.Experiment & key,'eye','prm','eyeClockTime','atTrialTime',0,'what','clocktime');
% Fit a line to translate eyelink time to nsTime. Even though some TRIALID
% events can be lost or delayed, a linear fit does a good job linking the
% two 
clockParms =  polyfit(trialIDTimeEyelink(~missingTrials),trialIDTimeNeurostim(~missingTrials),1); 
resid = polyval(clockParms,trialIDTimeEyelink(~missingTrials))-trialIDTimeNeurostim(~missingTrials);
slope = clockParms(1);
fprintf('Clock residuals: %3.3f ms +/- %3.3f ms, drift %3.3f ms/ms. %d missing TRIALID messages. \n',mean(resid),std(resid),slope-1,sum(missingTrials));
if (abs(slope-1)>0.5)
    error('Neurostim-Eyelink Clock skew larger than 0.5 ms/ms detected in edf file %s',filename);
end
if (mean(missingTrials)>0.1)
    warning('More than 10% TRIALIDs are missing from edf file %s. ',filename);
end
% Convert eyelink sample time to neurostim time
time  = polyval(clockParms,double(data.FSAMPLE.time))'/1000;
nrTimePoints= numel(time);
% Put the requested channels in signal
nrChannels =numel(parms.channel);
signal= nan(nrTimePoints,nrChannels);
for i=1:nrChannels
    signal(:,i) = data.FSAMPLE.(parms.channel{i})(parms.eye,:)';
end

%% Fill missing values 
MISSING = -32768; % From edf_data.h
signal(signal==MISSING)=NaN;
if isfield(parms,'fillmissing')
         fprintf('Filling missing %.1f%% of samples with %s ' ,100*mean(isnan(signal),'all'),parms.fillmissing{1});        
         signal = fillmissing(signal,parms.fillmissing{:});
end
%% Downsample
if isfield(parms,'downsample')
    R= ceil(parms.downsample/data.RECORDINGS(1).sample_rate);
    if R>1
        fprintf('Downsampling to %.0f Hz (decimate)...',parms.downsample);        
        tic
        nrTimePoints = ceil(nrTimePoints/R);
        tmp = nan(nrTimePoints,nrChannels);
        for ch = 1:nrChannels
            tmp(:,ch) =  decimate(signal(:,ch),R);
        end
        signal =tmp;
        time = linspace(time(1),time(end),nrTimePoints)';
        fprintf('Done in %d seconds.\n.',round(toc));
    end
end
% Information per channel
channelInfo =struct('name',parms.channel,'nr',num2cell(1:nrChannels));
% Overall recording info.
recordingInfo= data.RECORDINGS(end);
recordingInfo.header = data.HEADER;


end