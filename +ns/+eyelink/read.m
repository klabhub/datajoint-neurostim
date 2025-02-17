function[signal,time,channelInfo,recordingInfo] = read(key,parms)
% Read routine for ns.C to read from Eyelink edf data files, 
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
% .scale = Set to true to scale the raw pupil position measures to fractions of the
% camera image  (0,0) = center, (0.5,0.5) = top right. And pupil size
% between 0 and 1 with 1 the largest pupil that Eyelink captures.
%
% See Also fillmissing, decimate, ns.C
arguments
    key % The key of the Expriment table (File and CParm tuple)
    parms (1,1) struct  =struct% The preprocessing parameters
end

import ns.eyelink.*
% Fetch the file to read (ns.C has already checked that it exists)
filename = fullfile(folder(ns.Experiment &key),fetch1(ns.File &key,'filename'));


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
nrTrials =fetch1(ns.Experiment & key,'nrtrials');
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
time  = polyval(clockParms,double(data.FSAMPLE.time))';
nrSamples= numel(time);
% Put the requested channels in signal
nrChannels =numel(parms.channel);
signal= nan(nrSamples,nrChannels);
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
    R= ceil(data.RECORDINGS(1).sample_rate/parms.downsample);
    if R>1
        fprintf('Downsampling to %.0f Hz (decimate)...',parms.downsample);        
        tic
        nrSamples = ceil(nrSamples/R);
        tmp = nan(nrSamples,nrChannels);
        for ch = 1:nrChannels
            tmp(:,ch) =  decimate(signal(:,ch),R);
        end
        signal =tmp;
        time = linspace(time(1),time(end),nrSamples)';
        fprintf('Done in %d seconds.\n.',round(toc));
    end
end
if isfield(parms,'scale')
    fprintf('Scaling raw pupil coordinates\n');     
    % Raw pupil data (px py) are between -30e3 and + 30e3. We scale to this
    % to get position relative to the center (of the camera image), with
    % (0,0 at the center and (0.5,0.5) the top right corner.
    XYSCALE= 30e3;
    % The pupil area maxes out at 10000. Let's call that 1.
    PUPILSCALE = 10000;
    signal(:,ismember(parms.channel,{'px','py'})) = 0.5*signal(:,ismember(parms.channel,{'px','py'}))./XYSCALE;
    signal(:,ismember(parms.channel,{'pa'}))= signal(:,ismember(parms.channel,{'pa'}))./PUPILSCALE;
end
% Information per channel
channelInfo =struct('name',parms.channel,'nr',num2cell(1:nrChannels));
% Overall recording info.
recordingInfo= data.RECORDINGS(end);
recordingInfo.header = data.HEADER;
% Regular sampling so reduce time representation
time = double([time(1) time(end) nrSamples]);
% Reduce storage (ns.C.align converts back to double
signal  = single(signal);

%% Parse events to add as plugin parameter

startBlink = data.FEVENT([data.FEVENT.type]==3);
startBlinkTime = polyval(clockParms,double([startBlink.sttime]))';
trialStartTime = get(ns.Experiment & key,'cic','prm','firstFrame','what','clocktime');
nrBlinks = numel(startBlink);
startBlinkTrial = nan(nrBlinks,1);

for b= 1:nrBlinks 
    tmp = find(trialStartTime < startBlinkTime(b),1,'first');
    if isempty(tmp)
        startBlinkTrial(b)= 1;
    else
        startBlinkTrial(b) =tmp;
    end
end
startBlinkTrialTime = startBlinkTime- trialStartTime(startBlinkTrial);


plgTpl= fetch(ns.Experiment &key);
plgTpl.plugin_name ="edf";
 prmTpl  = plgTpl;   
 prmTpl.property_name = "startblink";
 if exists(ns.Plugin & plgTpl)
    delQuick(ns.PluginParameter &prmTpl);
    delQuick(ns.Plugin&plgTpl);
end


 prmTpl.property_value = [];
 prmTpl.property_time = startBlinkTrialTime; 
 prmTpl.property_nstime = startBlinkTime;
 prmTpl.property_trial = startBlinkTrial;
 prmTpl.property_type = 'Event';


 insert(ns.Plugin,plgTpl);
 insert(ns.PluginParameter,prmTpl)

end