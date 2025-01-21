function [signal,time] = CFilter(signal,time,parms)
% Generic downsampling and filtering function.
% Can be called from read functions (see intan.read for an example)
%
% parms.downsample - if defined, will downsample the signal to this
% frequency (in Hz),using decimate (and therefore an antialiasing filter). 
% Downsampling precedes filtering with parms.designfilt.
%
% parms.designfilt - defines filtering
% Filters are applied using filtfilt.
%
% Each member of the parms.designfilt struct contains one filter definition,
% all will be applied (successively, in order). 
%
% Example : 
% notchDesign = {'bandstopfir','FilterOrder',1000, ...
%                    'CutoffFrequency1',59,'CutoffFrequency2',61};
%
% parms = struct('designfilt',struct('notch',{notchDesign}));
[nrSamples,nrChannels] = size(signal);
sampleRate = 1./mode(diff(time));
   
%% Downsample first
if isfield(parms,'downsample')
    tic
    fprintf('Downsampling to %.0f Hz (decimate)...',parms.downsample);
    R=  round(sampleRate/parms.downsample);
    nrSamples = ceil(nrSamples/R);
    tmp = nan(nrSamples,nrChannels);
    for ch = 1:nrChannels
        tmp(:,ch) =  decimate(signal(:,ch),R);
    end
    signal =tmp;    
    time = linspace(time(1),time(end),nrSamples)';
    fprintf('Done in %d seconds.\n.',round(toc));
    sampleRate = parms.downsample;
end
%% Notch, Bandpass,etc.
if isfield(parms,'designfilt')
    fn = fieldnames(parms.designfilt);
    for i=1:numel(fn)
        tic;
        fprintf('Applying filter (designfilt.%s)...',fn{i})
        prms= parms.designfilt.(fn{i});
        d = designfilt(prms{:},'SampleRate',sampleRate);
        signal = filtfilt(d,signal);
        fprintf('Done in %d seconds.\n.',round(toc))
    end
end
