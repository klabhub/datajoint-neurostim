function [signal,time,channelInfo,recordingInfo] = spikeML(key,parms)
% Perform spike deconvolution on a Calcium fluorescence signal using the spikeML algorithm
% Deneux, T., Kaszas, A., Szalay, G., Katona, G., Lakner, T., Grinvald, A., Rozsa, B., & Vanzetta, I. (2016). Accurate spike estimation from noisy calcium signals for ultrafast three-dimensional imaging of large neuronal populations in vivo. Nature Communications 2016 7:1, 7(1), 1–17. https://doi.org/10.1038/ncomms12190
% Relies on the spikes and brick repositories: https://github.com/MLspike
%
% Wrapper by BK
% Sept 2024.
recordingInfo = '';
assert(exist('fn_structmerge','file'),"The brick repository must be on the path for spikeML");
assert(exist('spk_est.m','file'),"The spikes repository must be on the path for spikeML");

%% Check that the C table has a row with fluorescence for this key
% If that does not exist, the user needs to create it first, for instance
% by running suite2p preprocessing.
srcKey = key;
srcKey.ctag= "fluorescence";
srcC = ns.C & srcKey;
assert(exists(srcC),"SpikeML needs fluorescence data - generate these first in the ns.C table. For instance using suite2p (see sbx.read)");
time = fetch1(srcC,'time');

%% To store the drift estimates, create another row in ns.CParm
driftTag = key.ctag + "-drift";
if ~exists(ns.CParm & struct('ctag', driftTag))
    cPrm = fetch(ns.CParm & key,'*');
    cPrm.ctag =driftTag;
    cPrm.description = sprintf('This is the drift estimated during deconvolution for ctag = %s',key.ctag);
    insert(ns.CParm,cPrm);
end
driftKey = key;
driftKey.ctag = driftTag;


%% Fetch data and let spikeML do the work
[F,channel] = fetchn(ns.CChannel & srcKey,'signal','channel','LIMIT 1');
[t,dt] = sampleTime(srcC);
nrSamples = numel(t);
dt = dt/1000;% ms -> s
t = t /1000;
parms.parms.dt = dt;
nrChannels = numel(F);
signal = nan(nrSamples,nrChannels);
drift = cell(1,nrChannels);
channelInfo  =struct('nr',num2cell(channel));
for ch = 1:nrChannels
    if isfield(parms,'autocalibrate')
        % Do autocalibration
        pax = spk_autocalibration('par');
        pax.display = 'none'; % By default - can be overruled in parms.autocalibrate
        pax.dt = dt;
        fn = fieldnames(parms.autocalibrate);
        nrFields= numel(fn);
        for f=1:nrFields
            % Copy values from autocalibrate struct to pax
            % see spk_demo in mlspike/spikes
            % Values to set include
            % amin, amax, taumin,taumax, and saturation
            pax.(fn{f}) = parms.autocalibrate.(fn{f});
        end
        pax.mlspikepar.dographsummary = false;
        % perform auto-calibration
        fprintf('SpikeML: Autocalibrating Channel %d\n',ch)
        [tau,amp,sigma,events] = spk_autocalibration(F(ch),pax);
        calibratedParms = parms.parms; % Default from CParm
        % If events is empty, calibration was not possible, use defaults.
        if ~isempty(events)
            calibratedParms.tau = tau;
            calibratedParms.a = amp;
            calibratedParms.finetune.sigma= sigma;
        end
    else
        % Use parms as specified in CParm
        calibratedParms = parms.parms;
    end
    fprintf('SpikeML: Deconvolving Channel %d\n',ch)
    [spk,fit,drift{ch}] = spk_est(F{ch},calibratedParms);

    %% spk - spiketimes in s
    % Convert back to spike counts at the sample rate of the fluorescence
    signal(:,ch) = histcounts(spk,([t;t(end)+dt]-t(1)))';
    %% fit - fit of the F signal at each sample - used to estimate quality
    % and stored in channelInfo
    channelInfo(ch).quality = corr(fit,F{ch},Type="Pearson");
    channelInfo(ch).parms  = calibratedParms; % Store the parms that were used for this channel
    % If a previous populate call failed when inserting the deconvolved
    % spikes, the drift may still be in the database. Remove here, reinsert
    % the new values below.
    thisDriftTpl = mergestruct(driftKey,struct('channel',channel(ch)));
    if exists(ns.CChannel & thisDriftTpl)
        delQuick(ns.CChannel & thisDriftTpl);
    end
end


%% drift - drift of the F signal at each sample - add to separate ns.C row
% Create tuples and insert.
% Use a transaction because if something goes wrong here, none of the
% associated spiking data will be stored.
c= dj.conn;
c.startTransaction;
driftTpl = mergestruct(driftKey,...
    struct('time',time, ...
    'nrsamples',nrSamples,...
    'info',''));
insertIfNew(ns.C,driftTpl);

% Create tpl and insert
driftChannelTpl = mergestruct(driftKey,...
    struct('signal',num2cell(drift),...
    'channel',num2cell(channel),...
    'name','drift'));
insert(ns.CChannel,driftChannelTpl);
c.commitTransaction;

end