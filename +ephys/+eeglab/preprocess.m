function EEG  = preprocess(EEG,parms)
% Wrapper around EEGLAB preprocessing functions.
% Preprocess with EEGLab
% Fields of the parms.eeglab struct define the operations; they are
% executed in the sequence that they are defined in the parms.eeglab
% struct.
% 
% This called when populating hte ns.C table for an MFF file (from egi.read)
% Note that in egi.read the etc substruct of the EEG dataset is saved to keep track of 
% intermediate results outside the DJ database (And later inspection of
% results).
arguments
    EEG (1,1) struct
    parms (1,1) struct
end


fn = fieldnames(parms.eeglab);
for f= 1:numel(fn)
    switch fn{f}
        case 'ica'
            % Run ICA to identify components. The fields of the ica struct
            % are passed as parm/value pairs to pop_runica.           
            if ~isfield(parms.eeglab.ica,'pca')
                % The number of PCA components was not pre-specified.
                % Set it to the rank of the data. (Necessary for picard,
                % runica does this in pop_runica)
                nrSamplesToUse =min(3000,EEG.pnts); % Same as EEGLAB runica                
                X = double(EEG.data(:,1:nrSamplesToUse) - mean(EEG.data(:,1:nrSamplesToUse),2));                
                rnk =getrank(X);
                if rnk<size(X,1)
                    parms.eeglab.ica.pca = rnk;
                end
            end
            icaPV= namedargs2cell(parms.eeglab.ica);
            EEG = pop_runica(EEG,icaPV{:});
        case 'icaeog'
            % After running ICA, remove components based on EOG
            assert(~isempty(EEG.icaact),"Run ICA before removing ICA components");
            pv = setDefaults(parms.eeglab,'icaeog');
            % 1. Create hEOG by subtracting electrodes around the eye Left Outer Canthus (244) from Right Outer Canthus (237)
            hEOG = mean(EEG.data(pv.hEOG{1}, :),1,"omitmissing") - mean(EEG.data(pv.hEOG{2}, :),1,"omitmissing");
            vEOG = mean(EEG.data(pv.vEOG{1}, :),1,"omitmissing") - mean(EEG.data(pv.vEOG{2}, :),1,"omitmissing");
            % Calculate squared Pearson correlation coefficient
            v = abs(corr(EEG.icaact', vEOG'));
            h = abs(corr(EEG.icaact', hEOG'));
            r = max(v,h);
            % Find components where variance captured exceeds a threshold
            componentsToRemove = find(r> pv.rMin);
            fprintf('Components flagged for rejection: %s (r>=%.1f)\n', num2str(componentsToRemove'),min(r(componentsToRemove)));
            % Remove ICA components
            plotag = 0; keepcomp=0;
            EEG = pop_subcomp( EEG, componentsToRemove, plotag, keepcomp);
            % Add something to etc. for later retrieval
            EEG.etc.neurostim.icaeog = struct('r',r(componentsToRemove),'components',componentsToRemove);
        case 'icaeyetracker'
            assert(~isempty(EEG.icaact),"Run ICA before removing ICA components");
            pv = setDefaults(parms.eeglab,'icaeyetracker');

            preSamp  = round(pv.window  * EEG.srate / 1000);
            postSamp = round(pv.window* EEG.srate / 1000);
            nSamps   = preSamp + postSamp + 1;
            timesMs  = (-preSamp:postSamp) / EEG.srate * 1000;   % time axis (ms)

            EEG = ephys.eeglab.addEvents(EEG,pv.plugin,pv.events);
            assert(any(ismember(pv.events, {EEG.event.type})),"No %s events found in this dataset.",strjoin(pv.events));

            nEventTypes = numel(pv.events);
            % Variance explained by each ICA component (% of total EEG variance)

            nComp     = size(EEG.icaact, 1);
            varExplained = icaVarianceExplained(EEG);


            % Extract ETAs for each event type: [nComponents x nSamps x nEventTypes]
            eta = nan(nComp, nSamps, nEventTypes);
            nEventsByType = zeros(nEventTypes, 1);
            allEventTypes = erase(lower(string({EEG.event.type})), "'");
            for iEvt = 1:nEventTypes
                thisType = pv.events(iEvt);
                isThisType = allEventTypes == thisType;
                latencies = round([EEG.event(isThisType).latency]);
                ok = latencies > preSamp & latencies <= (EEG.pnts - postSamp);
                latencies = latencies(ok);
                nEventsByType(iEvt) = numel(latencies);
                if nEventsByType(iEvt) == 0
                    warning('No valid events found for type %s within epoch bounds.', thisType);
                    continue;
                end

                epochs = zeros(nComp, nSamps, nEventsByType(iEvt));
                for k = 1:nEventsByType(iEvt)
                    idx = latencies(k) - preSamp : latencies(k) + postSamp;
                    epochs(:, :, k) = EEG.icaact(:, idx);
                end
                eta(:, :, iEvt) = mean(epochs, 3);
            end

            % Baseline window [-200 -100] ms and response window [-100 200] ms
            baseWin = timesMs >= -pv.window & timesMs < -pv.window+pv.baseline;
            respWin = timesMs >= -pv.window+pv.baseline & timesMs <= pv.window;

            mu_base = squeeze(mean(eta(:, baseWin, :), 2));   % [nComp x nEventTypes]
            sd_base = squeeze(std(eta(:, baseWin, :), 0, 2)); % [nComp x nEventTypes]
            mu_resp = squeeze(mean(eta(:, respWin, :), 2));   % [nComp x nEventTypes]
            
            % z-score per component per event type
            zScore = (mu_resp - mu_base) ./ sd_base;          % [nComp x nEventTypes]

            % Remove component if any event crosses z-threshold and variance exceeds minimum.
            componentsToRemove = find(any(abs(zScore) > pv.minZ, 2) & (varExplained > pv.minVarPct))';

            % Add something to etc. for later retrieval
            eta_bc = eta - reshape(mu_base, [nComp 1 nEventTypes]);
            EEG.etc.neurostim.icaeyetracker = struct('eta',eta_bc,'z',zScore,'var',varExplained,'components',componentsToRemove);

            if pv.plot
                % Baseline-correct ETA with event-specific baseline               
                if isempty(componentsToRemove)
                    fprintf('No ICA components met removal criteria; skipping ETA/topoplot removal plots.\n');
                else
                    figure('Name','ETA of ICA components selected for removal');
                    clf
                    tiledlayout('flow');

                    nexttile;
                    hold on;
                    eventStyles =  {'-',':','-.'};
                    for iEvt = 1:nEventTypes
                        h = plot(timesMs, eta_bc(componentsToRemove, :, iEvt), 'LineWidth', 1.2,'LineStyle',eventStyles{iEvt});
                        for iLine = 1:numel(h)
                            h(iLine).DisplayName = sprintf('IC %d | %s', componentsToRemove(iLine), pv.events(iEvt));
                        end
                    end
                    h = xline(0, '--k');
                    h.DisplayName = 'Event Start';
                    xlabel('Time (ms)');
                    ylabel('Activation');
                    title 'Event Triggered Average (Baseline-Corrected)';

                    legend('Location', 'best');

                    % Plot 2D scalp topographies of removed ICA components.
                    for iComp = componentsToRemove
                        nexttile;
                        topoplot(EEG.icawinv(:, iComp), EEG.chanlocs(EEG.icachansind), ...
                            'electrodes', 'off', 'maplimits', 'absmax');
                        title(sprintf('IC %d | var=%.1f%% | max|z|=%.2f', iComp, varExplained(iComp), max(abs(zScore(iComp, :)), [], 'omitnan')));
                    end

                    saveas(gcf, fullfile(EEG.filepath, strrep(EEG.filename,'.mff','_ica.pdf')));
                    
                end
            end
            % Remove ICA components
            fprintf('%d components removed: %s (Var=%.1f%%)\n', numel(componentsToRemove),strjoin(string(componentsToRemove),'/'),sum(varExplained(componentsToRemove)));
            plotag = 0; keepcomp=0;
            EEG = pop_subcomp( EEG, componentsToRemove, plotag, keepcomp);
       
        case 'resample'
            if iscell(parms.eeglab.resample)
                % Passed verbatim to pop_resample
                % Must be {freq [Hz] ,fc,df}
                %   fc         - anti-aliasing filter cutoff (pi rad / sample)
                %                {default 0.9}
                %   df         - anti-aliasing filter transition band width (pi rad /
                %                sample) {default 0.2}
                % fc and df are optional
                resampleParms = parms.eeglab.resample;
            elseif isnumeric(parms.eeglab.resample) && isscalar(parms.eeglab.resample)
                % Only frequency specified.
                resampleParms = {parms.eeglab.resample};
            else
                error('parms.eeglab.resample must be either a scalar double (frequency) or a cell array with frequency,fc, and df. see pop_resample');
            end
            EEG = pop_resample(EEG, resampleParms{:});
        case 'zapline'
            if isstruct(parms.eeglab.zapline)
                if ~isfield(parms.eeglab.zapline,'noisefreqs')
                    % Oddly zapline does not look at line (50/60) noise
                    % by default. Force that default here
                    parms.eeglab.zapline.noisefreqs = 'line';
                end
                zapParms = namedargs2cell(parms.eeglab.zapline);
            elseif islogical(parms.eeglab.zapline) && parms.eeglab.zapline
                zapParms = {'noisefreqs','line'};
            end
            EEG = pop_zapline_plus(EEG, zapParms{:});
        case 'prep'
            % Use the PREP pipelin eegLab plugin for preprocessing
            % All unspecified values will be taken from the PREP defaults, except the
            % file paths (which we set to match the MFF file).
            if isstruct(parms.eeglab.prep)
                % The user specified parameters that differ from the PREP defaults
                % using a structure with structure fields.
            elseif islogical(parms.eeglab.prep) && parms.eeglab.prep
                % User had parms.prep =true
                % Use all PREP defaults
                %{
        % Not needed documentation only
        parms.prep.general.errorMsgs  = 'verbose';
        parms.prep.boundary.ignoreBoundaryEvents = false;
        parms.prep.resample.resampleOff = true;
        parms.prep.detrend.detrendChannels = 1:nrChannels;
        parms.prep.detrend.detrendType = 'high pass';
        parms.prep.detrend.detrendCutOff = 1;   % Hz
        parms.prep.detrend.detrendStepSize = 0.02; % Seconds
        %parms.prep.globaltrend =- not used.
        parms.prep.linenoise.lineNoiseMethod = 'clean';
        parms.prep.linenoise.lineNoiseChannels = 1:nrChannels;
        parms.prep.linenoise.Fs = EEG.srate;
        parms.prep.linenoise.lineFrequencies = 60:60:round(EEG.srate/2);
        parms.prep.linenoise.p  = 0.01;    
        parms.prep.linenoise.fScanBandwith = 2;
        parms.prep.linenoise.taperBandwith = 2;
        parms.prep.linenoise.taperWindowSize = 4;
        parms.prep.linenoise.taperWindowStep =1;
        parms.prep.linenoise.tau =100;
        parms.prep.linenoise.pad =0 ;
        parms.prep.linenoise.fPassBand = [0 EEG.srate/2];
        parms.prep.linenoise.maximumIterations = 10;
        parms.prep.reference.srate = EEG.srate;
        parms.prep.reference.samples = nrSamples;
        parms.prep.reference.robustDeviationThreshold = 5;
        parms.prep.reference.highFrequencyNoiseThreshold = 5;
        parms.prep.reference.correlationWindowSeconds =1;
        parms.prep.reference.correlationThreshold = 0.4;
        parms.prep.reference.badTimeThreshold = 0.01;
        parms.prep.reference.ransacOff = false;
        parms.prep.reference.ransacSampleSize = 50;
        parms.prep.reference.ransacChannelFraction = 0.25;
        parms.prep.reference.ransacCorrelationThreshold = 0.75;
        parms.prep.reference.ransacUnbrokenTime = 0.4;
        parms.prep.reference.ransacWindowSeconds = 5;
        parms.prep.reference.referenceType = 'robust';
        parms.prep.reference.interpolationOrder = 'post-reference';
        parms.prep.reference.meanEstimateType = 'median';
        parms.prep.reference.referenceChannels = 1:nrChannels;
        parms.prep.reference.evaluationChannels =1:nrChannels;
        % parms.prep.reference.channelLocations = EEG.chanlocs;    %Dont set these
        % parms.prep.reference.channelInformation= EEG.chaninfo;
        parms.prep.reference.maxReferenceIterations =4;
        parms.prep.reference.reportingLevel ='verbose';
        parms.prep.report.reportMode  = 'normal';
        parms.prep.report.summaryFilePath = './summary.pdf';
        parms.prep.report.sessionFilePath = '.report.html';
        parms.prep.report.consoleFID = 1;
        parms.prep.report.publishOn = true;
        parms.prep.report.errorMsgs = 'verbose';
        parms.prep.postprocess.keepFiltered = false;
        parms.prep.postprocess.removeInterpolatedChannels = false;
        parms.prep.postprocess.cleanupReference = false;
                %}
            end
            % Only change the reporting file paths if they have not been
            % specified already in the prep parms
            if ~isfield(parms.eeglab.prep,'report')
                parms.eeglab.prep.report = struct('reportMode','normal');
            end
            % Reporting does not work with an absolute file path due to some
            % weirdness in the prep pipeline
            % Temporarily cd
            here= pwd;
            cd(EEG.filepath)
            if ~isfield(parms.eeglab.prep.report,'summaryFilePath')
                parms.eeglab.prep.report.summaryFilePath  = ['.' filesep 'prep_summary.html'];
            end
            if ~isfield(parms.eeglab.prep.report,'sessionFilePath')
                parms.eeglab.prep.report.sessionFilePath  =  ['.' filesep strrep(EEG.filename,'.mff','_prep.pdf')];
            end
            try
                EEG = pop_prepPipeline(EEG, parms.eeglab.prep);
            catch me
                cd (here)
                rethrow(me)
            end
            % Keep only the channels that are not marked as stillNoisy
            EEG = pop_select(EEG, 'channel', setdiff(1:EEG.nbchan,EEG.etc.noiseDetection.stillNoisyChannelNumbers)');
            cd (here)
        case 'filt'
            EEG = pop_eegfiltnew(EEG, 'locutoff',parms.eeglab.filt.locutoff,'hicutoff',parms.eeglab.filt.hicutoff,'plotfreqz',0,'usefftfilt',true);
        otherwise
            error('Unknown eeglab preprocessing struct %s  \n',fn{f})
    end
    EEG = eeg_checkset(EEG);
end
end






function o = setDefaults(prms,fld)
arguments
    prms (1,1) struct
    fld (1,1) string
end
assert(isfield(prms,fld),'No %s field in the parms struct. Cannot continue.')

userParms= prms.(fld);
switch (fld)
    case 'icaeog'
        % Remove previously compute ICA components based on their
        % correlation with a synthetic EOG signal (constructed by
        % subtracting electrodes around the eye).
        defParms.hEOG = {237,244}; % Electrodes defining the horizontal EOG
        defParms.vEOG = {25,241}; % Vertical eog electrodes
        defParms.rMin = 0.5; % IC with a Pearson correlation with hEOG or vEOG that is bigger than this will be removed.
    case 'icaeyetracker'
        % Remove previously compute ICA components based on eye tracking
        % data
        defParms.window   = 200; % Use 200 ms before and after the event
        defParms.baseline = 100; % Use the first 100 ms of the window as the baseline
        defParms.events = "startblink";
        defParms.plugin = "edf";
        defParms.minZ       = 5;     % Remove components with an ETA that has a z-score higher than this
        defParms.minVarPct  = 1;  % Remove components only if they capture at least this amount of variance
        defParms.plot       = true;    % Plot ETA traces of components selected for removal


    otherwise
        error('No default settings for %s\n',fld);
end

userFn = string(fieldnames(userParms));
defFn = string(fieldnames(defParms));
extraneous = setdiff(userFn,defFn);
assert(isempty(extraneous),'There are extraneous fields in the %s struct (%s)',fld,strjoin(extraneous,'/'));

% Start with defaults, then overwrite the user-specified parms
o = defParms;
if ~isempty(userFn)
    for f=userFn(:)'
        o.(f) = userParms.(f);
        if iscellstr(o.(f)) || ischar(o.(f)) %#ok<ISCLSTR>
            o.(f) =  string(o.(f));
        end
    end
end
end

function tmprank2 = getrank(tmpdata)
% From EEGLAB - the cov/eig path seems necessary    
    tmprank = rank(tmpdata);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Here: alternate computation of the rank by Sven Hoffman
    %tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2)))); old code
    covarianceMatrix = cov(tmpdata', 1);
    [~, D] = eig (covarianceMatrix);
    rankTolerance = 1e-7;
    tmprank2=sum (diag (D) > rankTolerance);
    if tmprank ~= tmprank2
        tmprank2 = min(tmprank, tmprank2);
    end
end