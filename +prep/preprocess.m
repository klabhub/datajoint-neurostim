function [signal,time,result,channels] = preprocess(signal,time,parms,key)
% Generic downsampling and filtering function.
% Can be called from read functions (see intan.read for an example)
% signal = [nrSamples nrChannels]
% time   = [nrSamples 1] time in seconds
% channels = [1 nrChannels] - Channel numbers that were kept after
%                               preprocessing.
% The user passes a parms struct with any of the following fields.
% These filtering operations are applied in order (i.e. as ordered by
% fieldnames(parms). To apply the same filter multiple times, suffix the
% name with a number (1-9).
%
% parms.detrend  - a cell array of parameters passed to detrend.m
%                   Detrend the signal
% parms.decimate - a cell array of parameters passed to decimate
%                   Downsample the signal, using decimate (and therefore an antialiasing filter).
%                   If you prefer to specify a target frequency (instead of
%                   the downsamplng factor), use {'frequency',120} to
%                   decimate to 120Hz.
% parms.filtfilt -  each field is a cell array that defines a filter design (passed to
%                           designfilt), the filters are applied using
%                           filtfilt, sequentiallly in the order in which they are
%                           defined in parms.filtfilt
%
% parms.noisyChannels  - Find noise EEG channels using PREP criteria  (See
% find_noisy_channels.m)
%
% parms.reference  - Rereference all channels using the average, or
%                   laplacian.
%
% key - the key (a tuple) that this signal belongs too. Currently only used
% by the find_noisy_channels option.
%
% Example :
% First Reduce sampling by a factor of 4 (e.g. from 1kHz to 250Hz)
% parms.decimate = {4};
% Then filter out the 60 Hz noise followed by a low-pass
% parms.filtfilt.notch = {'bandstopfir','FilterOrder',1000, 'CutoffFrequency1',59,'CutoffFrequency2',61};
% parms.filtfilt.lowpass = {'lowpassfir','FilterOrder',2,'CutOffFrequency',120}
% and then removes a single linear trend
% parms.detrend  = {1,'Continuous', true};

if ismatrix(signal) % [samples channels]
    signal = permute(signal,[1 3 2]); % Add singleton trial dimension
end
[nrSamples,nrTrials,nrChannels] = size(signal);
channels = 1:nrChannels;  
if isfield(parms,"badChannel")
    badChannels   = parms.badChannel.channels;
    removeBadChannels = parms.badChannel.remove;
else
    badChannels   = [];
    removeBadChannels = false;
end                
samplingRate = 1./mode(diff(time));
fn =string(fieldnames(parms))';
result = struct('dummy',true);
for f=fn
    % In order to apply same transformations multiple times at different
    % orders the parms field can end with an identifier number, which is
    % removed at this step.
    match = regexp(f,'(?<command>\w+)(?<nr>\d+$)','names');
    if isempty(match)
        command =f;
    else
        command =match.command;
    end
    thisParms = parms.(f); % For ease of reference below
    if isfield(thisParms,'options')
        options =thisParms.options;
    else
        options = {};
    end
    tic
    % Developers: If a command changes not only signal, but also time,
    % or the badChannels list make sure to update those variables 
    % too , so that the  next step can continue with the correct information.
    switch command        
        case "baseline"
            isBaseline = do.ifinrange(time,thisParms.window);
            base = mean(signal(isBaseline,:,:), 1, "omitmissing");
            if isfield(thisParms,'decibel') && thisParms.decibel
                % Ratio as decibel
                signal = 10 * log10(signal./base);
            else
                % Subtract
                signal = signal - base;
            end
        case "resample"
            [signal,time] =resample(signal, thisParms.frequency, samplingRate,options{:});            
        case "zapline"
            assert(exist("clean_data_with_zapline_plus","file"),"Zapline requires an external toolbox. Get our fork at https://github.com/klabhub/zapline-plus.git and add it to your path")
            signal = reshape(signal,nrSamples*nrTrials,nrChannels);            
            signal = clean_data_with_zapline_plus(signal,round(samplingRate),thisParms);
            signal = reshape(signal,nrSamples,nrTrials,nrChannels);
        case "decimate"
            %% Downsampling using decimate
            % thisParms.frequency => specifies the target frequency
            targetRate = thisParms.frequency; % Set target frequency for decimation
            R = round(samplingRate/targetRate);
            fprintf('Downsampling from %.0f Hz to to %.0f Hz (decimate)...',samplingRate,targetRate);
            nrSamples = ceil(nrSamples/R);
            tmp = nan(nrSamples,nrTrials,nrChannels);
            for ch = 1:nrChannels
                for tr = 1:nrTrials
                    tmp(:,tr,ch) =  decimate(signal(:,tr,ch),R);
                end
            end
            signal =tmp;
            time = linspace(time(1),time(end),nrSamples)';           
        case "filtfilt"
            %% Notch, Bandpass,etc.
            % Any filter that can be designed with designfilt
            % and applied with filtfilt
            fprintf('Applying filter ...')
            d = designfilt(options{:},'SampleRate',samplingRate);
            signal = filtfilt(d,signal);
        case "detrend"
            %% Detrending using the detrend function
            fprintf('Detrending (%d order)...',thisParms.order)
            signal = detrend(signal,thisParms.order,"omitmissing",options{:});
        case "noisy_channels"
            fprintf('Finding noisy channels...\n')
            if isfield(parms, "layout")
                thisParms.ChannelLocations = parms.layout.ChannelLocations;
            end
            if isfield(thisParms, "epoch_buffer")
                epoch_mask = getTimepointsAroundTrials(ns.Experiment & key, time*1000, thisParms.epoch_buffer);
                thisParms = rmfield( thisParms,"epoch_buffer"); % No longer needed.
            else
                epoch_mask =  true(nrSamples*nrTrials,1);
            end
            thisParms.Fs = samplingRate;
            pvPairs = namedargs2cell(thisParms);
            tmpSignal = reshape(signal(epoch_mask,:,:),sum(epoch_mask)*nrTrials,nrChannels);
            result.noisy_channels = prep.find_noisy_channels(tmpSignal',pvPairs{:});
            result.noisy_channels.parameters = rmfield(result.noisy_channels.parameters, 'ChannelLocations'); % duplicate
            badChannels = union(badChannels,result.noisy_channels.all);            
        case "reference"
            fprintf('Re-referencing with %s mode',thisParms.method);
            signal = reshape(signal,nrSamples*nrTrials,nrChannels);
            if isempty(badChannels)
                referenceChannels = channels; % Use all.
            else
                if isfield(thisParms, "handleBads") && strcmpi(thisParms.handleBads, "interpolate")
                    %% Fix here, add option for interpParms{:} and chosing interp_func
                    signal = ns.interpolate_by_inverse_distance(signal, parms.layout.chanLocs, goodChannels, badChannels);
                    referenceChannels = channels; % After interp we use all channels
                else
                    referenceChannels = goodChannels; % Use only good channels (i.e. exclude bad)
                end
            end
            switch thisParms.method
                case "average"
                    reference = mean(signal(:,referenceChannels),2,"omitmissing");
                case "channel"
                    reference  = mean(signal(:,intersect(thisParms.channel,referenceChannels)),2,"omitmissing");
                case "laplacian"
                    reference = laplacian_reference(signal, parms.layout.neighbors,referenceChannels);
                otherwise
            end
            if any(all(isnan(reference),1))
                fprintf('All NaN reference signal... \n')
            end
            signal = signal-reference;
            signal = reshape(signal,nrSamples,nrTrials,nrChannels);
        case "badChannel"
            % Nothing to do - handled after everything is done
        otherwise
            % Not a defined preprocessing operation, just skip.
            fprintf('Skipping unknown preprocessing operation %s.\n',command);
    end
    fprintf('Done in %d seconds.\n',round(toc));

    % Update after each step - channels or time points may have changed
    [nrSamples,nrTrials,nrChannels] = size(signal);    
    goodChannels = setdiff(channels,badChannels);
    samplingRate = 1./mode(diff(time));
end

% If requested, remove the bad channels
if removeBadChannels
    signal(:,:,badChannels) = [];
    channels(badChannels) =[];
end
 
if nrTrials ==1
    % Back to original shape (remove singleton trial)
    signal = squeeze(signal);
end
end

%% Subtract mean of neighbors from each channel
function tmp = laplacian_reference(signal, neighbors,referenceChannels)
tmp = zeros(size(signal));
for i = 1:numel(neighbors)
    tmp(:, neighbors(i).channelnumber) =  mean(signal(:,intersect(neighbors(i).neighbors,referenceChannels)), 2,"omitmissing");
end
end
