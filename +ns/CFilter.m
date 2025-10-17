function [signal,time,result] = CFilter(signal,time,parms,expt)
% Generic downsampling and filtering function.
% Can be called from read functions (see intan.read for an example)
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
% expt - the experiment that this signal belongs too. Currently only used
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

[nrSamples,nrChannels] = size(signal);
sampleRate = 1./mode(diff(time));
if isfield(parms,'badElectrodes') 
    badElectrodes = parms.badElectrodes;
else
    badElectrodes= [];
end
fn =string(fieldnames(parms))';
result = struct([]);
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
    tic            
    switch command
        case "zapline"
            assert(exist("clean_data_with_zapline_plus","file"),"Zapline line noise removal requires the zapline-plus toolbox. \n Get it on github and add it to your path")            
            signal = clean_data_with_zapline_plus(signal,round(sampleRate),thisParms);
        case "decimate"
            %% Downsampling using decimate
            % thisParms.frequency => specifies the target frequency   
            targetRate = thisParms.frequency; % Set target frequency for decimation
            R = round(sampleRate/targetRate);                            
            fprintf('Downsampling from %.0f Hz to to %.0f Hz (decimate)...',sampleRate,targetRate);                        
            nrSamples = ceil(nrSamples/R);
            tmp = nan(nrSamples,nrChannels);
            for ch = 1:nrChannels
                tmp(:,ch) =  decimate(signal(:,ch),R);
            end
            signal =tmp;
            time = linspace(time(1),time(end),nrSamples)';            
        case "filtfilt"
            %% Notch, Bandpass,etc. 
            % Any filter that can be designed with designfilt
            % and applied with filtfilt
            ffFn = fieldnames(thisParms);
            for i=1:numel(ffFn)
                tic;
                fprintf('Applying filter (designfilt.%s)...',ffFn{i})
                prms= thisParms.(ffFn{i});
                d = designfilt(prms{:},'SampleRate',sampleRate);
                signal = filtfilt(d,signal);                
            end
        case "detrend"
            %% Detrending using the detrend function                
            fprintf('Detrending (%d)...',thisParms{1})
            signal = detrend(signal,thisParms{:});
           
        case "noisy_channels"         
            fprintf('Finding noisy channels...\n')
            if isfield(parms, "layout")
                thisParms.ChannelLocations = parms.layout.ChannelLocations;
            end                               

            if isfield(thisParms, "epoch_buffer")
                    epoch_mask = ns.getTimepointsAroundTrials(expt, time*1000, thisParms.epoch_buffer);
                    thisParms = rmfield( thisParms,"epoch_buffer"); % No longer needed.
            else
                    epoch_mask =  ones(1,size(signal,1))==1;
            end
            thisParms.Fs = sampleRate; 
            pvPairs = namedargs2cell(thisParms);
            result.noisy_channels = ns.prep.find_noisy_channels(signal(epoch_mask,:)',pvPairs{:});
            result.noisy_channels.parameters = rmfield(result.noisy_channels.parameters, 'ChannelLocations'); % duplicate
            badElectrodes = union(badElectrodes,result.noisy_channels.all);            
        case "reference"
            fprintf('Re-referencing with %s mode',thisParms.method);            
            if ~isempty(badElectrodes) 
                if isfield(thisParms, "handleBads") && strcmpi(thisParms.handleBads, "interpolate")
                    %% Fix here, add option for interpParms{:} and chosing interp_func
                    signal = ns.interpolate_by_inverse_distance(signal, parms.layout.chanLocs, ...
                                                setdiff(1:nrChannels, badElectrodes), badElectrodes);
                    referenceChannels = 1:nrChannels;
                else 
                    referenceChannels = setdiff(1:nrChannels,badElectrodes); 
                end
            else
                referenceChannels = 1:nrChannels;
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
        otherwise
            % Not a defined filter operation, just skip.
    end
    fprintf('Done in %d seconds.\n',round(toc));

    % Update after each step (nrChannels does not change, but nrSamples
    % can)
    [nrSamples,nrChannels] = size(signal);
    sampleRate = 1./mode(diff(time));   
end
end

%% Subtract mean of neighbors from each channel
function tmp = laplacian_reference(signal, neighbors,referenceChannels)
    tmp = zeros(size(signal));
    for i = 1:numel(neighbors)            
        tmp(:, neighbors(i).channelnumber) =  mean(signal(:,intersect(neighbors(i).neighbors,referenceChannels)), 2,"omitmissing");    
    end
end
