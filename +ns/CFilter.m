function varargout = CFilter(signal,time,parms)
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
hasNoisyChannels = isfield(parms,'badElectrodes') && ~isempty(parms.badElectrodes);
varargout = cell(1, nargout);
fn =string(fieldnames(parms))';
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

    switch command
        case "decimate"
            %% Downsampling using decimate
            tic
            if iscell(parms.(f)) && strcmpi(parms.(f){1},"frequency")
                targetRate = parms.(f){2};
                R = round(sampleRate/targetRate);                
            else
                R=  parms.(f); % parms.decimate is the R factor
                targetRate =  sampleRate/R;
            end
            fprintf('Downsampling from %.0f Hz to to %.0f Hz (decimate)...',sampleRate,targetRate);                        
            nrSamples = ceil(nrSamples/R);
            tmp = nan(nrSamples,nrChannels);
            for ch = 1:nrChannels
                tmp(:,ch) =  decimate(signal(:,ch),R);
            end
            signal =tmp;
            time = linspace(time(1),time(end),nrSamples)';
            fprintf('Done in %d seconds.\n',round(toc));            
        case "filtfilt"
            %% Notch, Bandpass,etc. 
            % Any filter that can be designed with designfilt
            % and applied with filtfilt
            fn = fieldnames(parms.(f));
            for i=1:numel(fn)
                tic;
                fprintf('Applying filter (designfilt.%s)...',fn{i})
                prms= parms.(f).(fn{i});
                d = designfilt(prms{:},'SampleRate',sampleRate);
                signal = filtfilt(d,signal);
                fprintf('Done in %d seconds.\n',round(toc))
            end
        case "detrend"
            %% Detrending using the detrend function    
            tic;
            fprintf('Detrending (%d)...',parms.(f){1})
            signal = detrend(signal,parms.(f){:});
            fprintf('Done in %d seconds.\n',round(toc));

        case "noisy_channels"

            tic;
            fprintf('Finding noisy channels...\n')
            parmsN = parms.(f);

            if isfield(parms, "layout")
                parmsN.ChannelLocations = parms.layout.ChannelLocations;
            end
            
            ep_mask = parms.(f).epoch_mask;
            parmsN = rmfield(parmsN,"epoch_mask");
            if isfield(parmsN,"epoch_buffer")
                parmsN = rmfield(parmsN, "epoch_buffer");
            end
            parmsN = namedargs2cell(parmsN);
            varargout{3} = ns.prep.find_noisy_channels(signal(ep_mask,:)', parmsN{:});
            varargout{3}.parameters = rmfield(varargout{3}.parameters, 'ChannelLocations'); % duplicate
            hasNoisyChannels = ~isempty(varargout{3}.all);
            parms.badElectrodes = unique(horzcat(parms.badElectrodes, varargout{3}.all));
            fprintf('Done in %d seconds.\n',round(toc));

        case "reference"
            ref_opts = parms.(f);
            isBadsArg = ismember(ref_opts, "bads");
            if any(isBadsArg) && hasNoisyChannels
                handle_bads = ref_opts{find(isBadsArg) + 1};                
            elseif hasNoisyChannels
                handle_bads = "exclude";
            else
                handle_bads = [];
            end

            if ~isempty(handle_bads) && strcmp(handle_bads, "interpolate")
                %% Fix here, add option for interpParms{:} and chosing interp_func
                signal = ns.interpolate_by_inverse_distance(signal, parms.layout.chanLocs, ...
                    setdiff(1:size(signal,1), parms.badElectrodes), parms.badElectrodes);
            end

            switch parms.(f){1}
                case "average"
                    disp('Applying rereferencing to the average channel activity.');
                    signal = signal - mean(signal,2, 'omitnan', true);
                case "channel"
                case "laplacian"
                    disp('Applying Laplacian re-referencing channel activity.');
                    signal = laplacian_reference(signal, parms.layout.neighbors);                   
                otherwise
                    continue
            end
        otherwise
            % Not a defined filter operation, just skip.
    end
    % Update after the filter step
    [nrSamples,nrChannels] = size(signal);
    sampleRate = 1./mode(diff(time));
    varargout{1} = signal;
    varargout{2} = time;
end
end

%% Subtract mean of neighbors from each channel
function tmp = laplacian_reference(signal, neighbors)
    tmp = zeros(size(signal));
    for iCh = 1:width(neighbors)    
        tmp(:, iCh) = signal(:,iCh) - mean(signal(:,neighbors(iCh).neighbors), 2);    
    end
end
