function EEG  = preprocess(EEG,parms)
% Wrapper around EEGLAB preprocessing functions. 
% Preprocess with EEGLab
% Fields of the parms.eeglab struct define the operations; they are
% executed in the sequence that they are defined in the parms.eeglab
% struct.
arguments
    EEG (1,1) struct
    parms (1,1) struct
end


fn = fieldnames(parms.eeglab);
for f= 1:numel(fn)
    switch fn{f}
        case 'ica'

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
               %parms.eeglab.prep =struct('report',struct);
            end
            % Only change the reporting file paths if they have not been
            % specified already in the prep parms
            fldr = folder(ns.Session & key);
            if ~isfield(parms.eeglab.prep,'report')
                parms.eeglab.prep.report = struct('reportMode','normal');
            end
            % Reporting does not work with an absolute file path due to some
            % weirdness in the prep pipeline
            % Temporarily cd
            here= pwd;
            cd(fldr)
            if ~isfield(parms.eeglab.prep.report,'summaryFilePath')
                parms.eeglab.prep.report.summaryFilePath  = ['.' filesep 'prep_summary.html'];
            end
            if ~isfield(parms.eeglab.prep.report,'sessionFilePath')
                parms.eeglab.prep.report.sessionFilePath  =  ['.' filesep char(strrep(key.filename,'.mff','_prep.pdf'))];
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
end
end