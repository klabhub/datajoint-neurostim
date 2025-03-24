function [spk,time,channelInfo,recordingInfo] = spikeML(key,parms)
% Perform spike deconvolution on a Calcium fluorescence signal using the spikeML algorithm
% Deneux, T., Kaszas, A., Szalay, G., Katona, G., Lakner, T., Grinvald, A., Rozsa, B., & Vanzetta, I. (2016). Accurate spike estimation from noisy calcium signals for ultrafast three-dimensional imaging of large neuronal populations in vivo. Nature Communications 2016 7:1, 7(1), 1–17. https://doi.org/10.1038/ncomms12190
% Relies on the spikes and brick repositories: https://github.com/MLspike
%
% Wrapper by BK
%
% parms contains instructions on the source rows in C (i.e., the
% fluorescence and the neuropil background), and the number of parallel
% workers:
%
% parms.nrWorkers = 5;  % parfor loop over channels in spikeML
% parms.ctagF = "fluorescence"; % The ctag in ns.C that has the fluorescence data (should already be populated)
% parms.ctagBg = "neuroopil"; % The ctag in ns.C that has the fluorescence data (should already be populated)
%% The .deconv substruct defines deconvolution parameters.
% Note that these parameters follow the definition of tps_mlspikes
% a , tau, and sigma can be determined by autocalibration, but if that fails
% (no isolated events) then these defaults weill be used:
% parms.deconv = tps_mlspikes('par'); % Initialize parameters for the MLspike algorithm.
% parms.deconv.a = 0.07; % Expected amplitude.
% parms.deconv.tau = 1.3; % Gcamp6s is 1.3s
% parms.deconv.algo.nspikemax = 4; % Set the maximum number of spikes expected within one time window.
% parms.deconv.pnonlin = [0.85 -0.006];  % Deneux values for Gcamp6s.
% parms.deconv.drift.parameter = 0.01; % Drift parameters.
% parms.deconv.drift.baselinestart = 1; % Specify the starting frame for baseline correction.
% parms.deconv.dographsummary = false; % Disable graphical summary of the spike inference.
% parms.deconv.display = 'no'; % Disable display outputs during spike inference.
%%  Define autocalibration parameters.
% Autocalibration (if defined, by creating the .autocalibrate field in the parms struct)
% can (if successful) overrule the a, tau, and sigma parameters for the deconvolution.
% Here too, the definitions come from the MLspike toolbox:
% parms.autocalibrate = spk_autocalibration('par'); % Initialize parameters for autocalibration.
% parms.autocalibrate.amin = 0.03; % Set minimum amplitude for spike amplitude estimation.
% parms.autocalibrate.amax = 0.19; % Set maximum amplitude for spike amplitude estimation.
% parms.autocalibrate.taumin = 1.52; % Set minimum decay time constant for the calcium indicator.
% parms.autocalibrate.taumax = 2.22; % Set maximum decay time constant for the calcium indicator.
% parms.autocalibrate.driftparam = 0.01; % Set the drift parameters, possibly to account for baseline drift in fluorescence signal.
% parms.autocalibrate.pnonlin = [0.85 -0.006];  % Deneux values for Gcamp6s -  Set parameters for non-linear processing, if applicable.
% parms.autocalibrate.mlspikepar.dographsummary = false; % Disable graphical summary of the autocalibration process.
% parms.autocalibrate.display = 'no'; % Disable display outputs during autocalibration.
%%
% Having defined the instructions (parms) we can add it as a CParm
% spikeMLPrep = struct('ctag','spikesMLAutoCal',...       % Name of this C
%                   'description','Deconvolution using spike ML with autocalibration',...
%                    'extension','.sbx',...  % Files to process
%                    'fun','sbx.spikeML',...  % function to call.
%                    'parms',parms);   % parameters to pass
% insertIfNew(ns.CParm, spikeMLPrep);
% Then:
% populate(ns.C,'ctag="spikesMLAutoCal"') wil populate the ns.C table with spikes
% for all sbx files. Files that don't have a fluorescence or neuropil entry in ns.C will fail.
%
% NOTE
%  The newly created deconvolved spike rows in ns.C are not linked to their "parent"
%  fluorescence and will not be deleted/modified  if their parents are deleted/modified.
%
% Sept 2024.


recordingInfo = '';
assert(exist('fn_structmerge','file'),"The brick repository must be on the path for spikeML");
assert(exist('spk_est.m','file'),"The spikes repository must be on the path for spikeML");
assertRequiredMembers(parms,["ctagF" "ctagBg" "deconv"],'sxb.spikeML');

%% Check that the C table has a row with fluorescence for this key
% If that does not exist, the user needs to create it first, for instance
% by running suite2p preprocessing.
fKey = key;
fKey.ctag = parms.ctagF; % Identifies the source of the fluorescence data by its ctag
fSrc = ns.C & fKey;
assert(exists(fSrc),"SpikeML needs %s data - generate these first in the ns.C table. For instance using suite2p (see sbx.read)",fKey.ctag);

bgKey = key;
bgKey.ctag = parms.ctagBg; % Identifies the source of the fluorescence data by its ctag
bgSrc = ns.C & bgKey;
assert(exists(bgSrc),"SpikeML needs %s data - generate these first in the ns.C table. For instance using suite2p (see sbx.read)",bgKey.ctag);

time = fetch1(fSrc,'time');

%% To store the drift estimates, create another row in ns.CParm
% but we exclude all files so that it will never lead to a
% maketuple call from populate.
driftTag = key.ctag + "-drift";
if ~exists(ns.CParm & struct('ctag', driftTag))
    cPrm = fetch(ns.CParm & key,'*'); % Copy the CParm from the MLSpike entry
    cPrm.ctag =driftTag; % Rename the ctag
    cPrm.exclude = "%"; % Exclude every file
    cPrm.description = sprintf('This is the drift estimated during deconvolution for ctag = %s',key.ctag);
    insert(ns.CParm,cPrm);
end
driftKey = key;
driftKey.ctag = driftTag;


%% Setup data structs to store results 
% 
[t,dt] = sampleTime(fSrc); %(should match neuropil timing)
nrSamples = numel(t);
dt = dt/1000;% ms -> s
t = t /1000;
parms.deconv.dt = dt;
if isfield(parms,'restrictRoi')
    % ROI restriction - allow user to specify something like
    % 'pcell>0.1 AND radius> 2.5'
    restrictRoi = proj(sbx.PreprocessedRoi & parms.restrictRoi,'roi->channel');
else
    restrictRoi  = true; % No restrictons
end
fChannels= fSrc*(ns.CChannel & restrictRoi);
bgChannels = bgSrc*(ns.CChannel & restrictRoi);
nrChannels = count(fChannels);

spk = nan(nrSamples,nrChannels);
drift = cell(1,nrChannels);
channelInfo  =repmat(struct('nr',nan,'nanFrac',nan,'quality',nan,'autoCal',nan,'parms',struct),[nrChannels 1]);

%Start parpool if requested
if isfield(parms,'nrWorkers') && parms.nrWorkers>0
    pool = gcp("nocreate");
    if isempty(pool)
        pool = parpool(parms.nrWorkers);
    else 
        parms.nrWorkers = max(parms.nrWorkers,pool.NumWorkers);        
    end
else
    pool = [];
end
fprintf('Deconvolving %d channels \n',nrChannels)
%% Loop for/parfor per channel
if isempty(pool)
    for  ch = 1:nrChannels
        [spk(:,ch),drift{ch}, channelInfo(ch)] = loopBody(fChannels,bgChannels,parms,driftKey,t,dt,ch);
    end
else
    parfor  (ch = 1:nrChannels,parms.nrWorkers)
        dj.conn; % Need to refresh connection in each worker
        [spk(:,ch),drift{ch}, channelInfo(ch)]  = loopBody(fChannels,bgChannels,parms,driftKey,t,dt,ch);
    end
end

%% Create tuples and insert.
% Drift
driftTpl = mergestruct(driftKey,...
    struct('time',time, ...
    'nrsamples',nrSamples,...
    'info',''));
insertIfNew(ns.C,driftTpl);
% Drift per channel
driftChannelTpl = mergestruct(driftKey,...
    struct('signal',drift',...
    'channel',{channelInfo.nr}',...
    'name',"drift"));
insert(ns.CChannel,driftChannelTpl);

% The spikes are returened to the calling function (ns.C/makeTuples) and
% inserted there.
end

%% Function that does the deconvolution for one channel 
% This is called from eitehr a for (nrWorkers=0) or parfor loop
function [spk,drift, channelInfo] = loopBody(fChannels,bgChannels,parms,driftKey,t,dt,ch)

[F,channel] = fetchn(fChannels,'signal','channel',['ORDER BY channel LIMIT 1 OFFSET ' num2str(ch-1)]);
[Bg] = fetchn(bgChannels,'signal',['ORDER BY channel LIMIT 1 OFFSET ' num2str(ch-1)]);
signal = (F{1}-Bg{1}); % Subtract the neuropil background from the fluorescence.
isNaN = isnan(signal);
signal(isNaN) = 0; % Could remove samples instead or linearly interpolate, but this should be rare (missing F)
if isfield(parms,'autocalibration') && ~isempty(fieldnames(parms.autocalibration))
    % Do autocalibration
    pax = spk_autocalibration('par'); % Get defaults, then overrule with parms.autocalibrate
    pax.display = 'none';
    pax.mlspikepar.dographsummary = false;
    fn = fieldnames(parms.autocalibration);
    nrFields= numel(fn);
    for f=1:nrFields
        % Copy values from autocalibrate struct to pax
        % see spk_demo in mlspike/spikes
        % Values to set include
        % amin, amax, taumin,taumax, and saturation
        pax.(fn{f}) = parms.autocalibration.(fn{f});
    end
    pax.dt = dt; % Dont allow overrule by autocalibrate.
    % perform auto-calibration
    fprintf('SpikeML: Autocalibrating Channel %d\n',channel)
    [tau,amp,sigma,events] = spk_autocalibration(signal,pax);
    calibratedParms = parms.deconv; % Default from CParm
    % If events is empty, calibration was not possible, use defaults.
    if isempty(events)
        autoCal = false;
    else
        calibratedParms.tau = tau;
        calibratedParms.a = amp;
        calibratedParms.finetune.sigma= sigma;
        autoCal = true;
    end
else
    % No autocalibration : use parms as specified in CParm
    calibratedParms = parms.deconv;
    autoCal = false;
end
fprintf('SpikeML: Deconvolving Channel %d\n',channel)
[thisSpk,fit,drift] = spk_est(signal,calibratedParms);

%% spk - spiketimes in s
% Convert back to spike counts at the sample rate of the fluorescence
spk = histcounts(thisSpk,([t;t(end)+dt]-t(1)))';
%% fit - fit of the F signal at each sample - used to estimate quality
% and stored in channelInfo
channelInfo.nr = channel;
channelInfo.quality = corr(fit(~isNaN),signal(~isNaN),Type="Pearson");
channelInfo.nanFrac = mean(isNaN);
channelInfo.parms  = calibratedParms; % Store the parms that were used for this channel
channelInfo.autoCal   = autoCal;
% If a previous populate call failed when inserting the deconvolved
% spikes, the drift may still be in the database. Remove here, reinsert
% the new values below.
thisDriftTpl = mergestruct(driftKey,struct('channel',channel));
if exists(ns.CChannel & thisDriftTpl)
    delQuick(ns.CChannel & thisDriftTpl);
end
end