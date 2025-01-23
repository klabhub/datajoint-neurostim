function [signal,time] = CFilter(signal,time,parms)
% Generic downsampling and filtering function.
% Can be called from read functions (see intan.read for an example)
% The user passes a parms struct with any of the following fields.
% These filtering operations are applied in order (i.e. as ordered by
% fieldnames(parms)
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

fn =string(fieldnames(parms))';
for f=fn
    switch f
        case "decimate"
            %% Downsampling using decimate
            tic
            if numel(parms.decimate)==2 & strcmpi(parms.decimate{1},"frequency")
                targetRate = parms.decimate{2};
                R = round(sampleRate/targetRate);                
            else
                R=  parms.decimate{1}; % First input to decimate is the R factor
                targetRate =  sampleRate/R;
            end
            fprintf('Downsampling to %.0f Hz (decimate)...',targetRate);                        
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
            fn = fieldnames(parms.filtfilt);
            for i=1:numel(fn)
                tic;
                fprintf('Applying filter (designfilt.%s)...',fn{i})
                prms= parms.filtfilt.(fn{i});
                d = designfilt(prms{:},'SampleRate',sampleRate);
                signal = filtfilt(d,signal);
                fprintf('Done in %d seconds.\n',round(toc))
            end
        case "detrend"
            %% Detrending using the detrend function    
            tic;
            fprintf('Detrending (%d)...',parms.detrend{1})
            signal = detrend(signal,parms.detrend{:});
            fprintf('Done in %d seconds.\n',round(toc))
        otherwise
            % Not a defined filter operation, just skip.
    end
    % Update after the filter step
    [nrSamples,nrChannels] = size(signal);
    sampleRate = 1./mode(diff(time));
end
end
