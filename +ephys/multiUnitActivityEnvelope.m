function [signal,muaeTime]=multiUnitActivityEnvelope(signal,time,pv)
% function [muae,muaeTime]=ana2muae(muae,time,...)
%
% Implementation of Super & Roelfsema 2005
% http://www.ncbi.nlm.nih.gov/pubmed/15581712
% In short:
% Bandpass filter the analog channel signal at 750-5000 Hz, rectify it and
% finally low pass filters at 500 Hz, resulting in the so called Multi Unit
% Activity Envelope (MUAE).
%
% Different values for low and high pass filtering, and the type of
% filtering, are controlled with the parms input argument
%
% INPUT:
%   signal:
%       analog signal, numeric. If this is a matrix, the conversion is done
%       for each columns. Must be double. NaNs will be set to zero before
%       filtering.
%   time:
%       time in seconds at which the signal was obtained, the sampling rate of the input
%       signal is determined from this. Must be a double.
% Optional Parameter/Value pairs
%   'parms':
%       struct with filter settings. Its format is shown here (which are
%       also the default settings).
%           struct('filt1',struct('lopHz',5000,'hipHz',750,'type','butterbandpass','order',4), ...
%                  'filt2',struct('lopHz',500,'hipHz',[],'type','butterlowpass','order',4), ...
%                   'targetHz',1000
%                   'maxMemory',4
%                    'useGPU',false); The maxmimum memory (in GB) available for the filtering operation.
%                                   Filtering will be done in chunks of this size. Default is 4 GB.
%
% OUTPUT:
%   muae:
%       multi unit activity envelope,
%   tmuae:
%       time axis of the muae
%
% Jacob, 2011-12-05
% BK Revised for DJ June 2023

% The input is called muae like the output to not
% keep two copies of this potentially large array in RAM.
arguments
    signal (:,:) double      % The raw data
    time (:,1) double    % The time in seconds at which the samples were obtained
    pv.parms (1,1) struct = struct('filt1',struct('lopHz',5000,'hipHz',750,'type','butterbandpass','order',4), ...
        'filt2',struct('lopHz',500,'hipHz',[],'type','butterlowpass','order',4), ...
        'targetHz',1000,...
        'maxMemory', 4,...
        'useGPU',false) % The memory available to do the filtering. In GB
end


if isrow(signal)
    makeColumn = true;
    signal=signal';
else
    makeColumn = false;
end

% Assuming that the
samplingRate =unique(round(1./diff(time),2));
assert(isscalar(samplingRate),"This sampling rate has more than 2 decimals??")

nannedData = sparse(isnan(signal));
signal(nannedData) = 0;

halfsampleHz=samplingRate/2;
if mod(pv.parms.filt1.order,2) || mod(pv.parms.filt2.order,2)
    error('filter orders must be even because filtfilt is used for filtering');
else
    halforder1=pv.parms.filt1.order/2;
    halforder2=pv.parms.filt2.order/2;
end

% Settings to do filtering piecemeal:
[nrSamples,nrSignals]=size(signal);
if pv.parms.maxMemory==0
    chunkSamples = nrSamples;
else
    chunkSamples = min(nrSamples,floor((pv.parms.maxMemory*1e9/(nrSamples*nrSignals*8))*nrSamples));
end
marginSamples = samplingRate; % One second margin overlap should be enough

% FILTER 1
% First, make a filter using the requested method (e.g., butter), then
% apply it using filtfilt, which, unlike filter, maintains phase. Because
% filtfilt filters twice (in two directions) we need to half the order of
% the filter. (On 2011-12-07 I checked that filter(b,a,sig) is less strong
% than filtfilt(b,a,sig) and found this to be true.)

switch pv.parms.filt1.type
    case 'butterbandpass'
        [b,a]=butter(halforder1/2,[pv.parms.filt1.hipHz pv.parms.filt1.lopHz]/halfsampleHz); % not only filtfilt, but also bandpass butter doubles the order
    case 'butterlowpass'
        [b,a]=butter(halforder1,[pv.parms.filt1.lopHz]/halfsampleHz,'low');
    case 'butterhighpass'
        [b,a]=butter(halforder1,[pv.parms.filt1.hipHz]/halfsampleHz,'high');
    otherwise
        error(['Unknown filtertype: ' pv.parms.filt1.type ]);
end
filtfiltpiecemeal;

% RECTIFICATION
signal=abs(signal);

% FILTER 2
switch pv.parms.filt2.type
    case 'butterbandpass'
        [b,a]=butter(halforder2,[pv.parms.filt2.hipHz pv.parms.filt1.lopHz]/halfsampleHz);
    case 'butterlowpass'
        [b,a]=butter(halforder2,[pv.parms.filt2.lopHz]/halfsampleHz,'low');
    case 'butterhighpass'
        [b,a]=butter(halforder2,[pv.parms.filt2.hipHz]/halfsampleHz,'high');
    otherwise
        error(['Unknown filtertype: ' pv.parms.filt2.type ]);
end
filtfiltpiecemeal;

% DOWNSAMPLE
skipstep=samplingRate/pv.parms.targetHz;
if mod(skipstep,1)~=0
    error('samplingRate must be a multiple of targetHz. If this is unacceptable, look into using resample instead of downsample.');
    % If this turns out to be a
    % problem, the code could be adapted to use the function resample.m to
    % reduce the maue length. However, because the result of this has a length
    % of ceil(size(ana,1)*(convset.filt2.lopass*2/samplingRate). In other
    % words, because of rounding error the realized samplerate of muae may not
    % exactly be convset.filt2.lopass*2. Also, it is not as straight forward
    % to make the corresponding time axis. Most importantly, however, the
    % function resample perfoms a low-pass filtering of it's own, making the
    % process more complicated to describe accurately.

    % use this code for domwsampling if the inputlength/outputlength ratio is
    % not a natural number, and write something to make a tmuae
    %muae=resample(TMP,convset.muae.targetsampleHz,samplingRate); % reduce
    %number of samples by factor convset.muae.targetsampleHz/samplingRate
    %muaeHz=length(muae)/length(ana)*samplingRate; % because of rounding errors,
    %use this as the realized muae sampling hz
    % note resample also applies a low pass filter to prevent aliasing, this
    % may make the apply filter 2 step unnecessary
end
signal=downsample(signal,skipstep);
muaeTime=downsample(time,skipstep);
nannedData = downsample(nannedData,skipstep);
signal(nannedData) = NaN; % Put them back


%% Nesting the piecemeal filter function to reduce memory requirements

    function filtfiltpiecemeal
        % Wrapper around filtfilt that does the filtering piecemeal to reduce
        % peak-RAM usage.
        %
        % INPUT:
        %   b,a,X: as in filtfilt
        %   chunksamples: size of chunks that signal vector array X will be split
        %       in (e.g. 3600/4*samplerateHz to split in 15 minutes intervals
        %   marginsamples: size in samples of margins around chunks to remove
        %       edge effects. I found that over 4 hour data and a lowpass of 600, a
        %       margin of a second results in identical results to filtfilt.
        %       Perhaps with lower low-pass frequencies a larger margin is needed.
        %
        % OUTPUT:
        %   X: filtered output of the exact length as input X
        %
        % With proper values for marginsamples, the output is exactly the same as
        % filtfilt.
        %
        % EXAMPLE
        %   disp('Open taskmgr/performance to monitor memory usage over time')
        %   samplinghz=25000;
        %   nsamp=samplinghz*3600*4; %  4 hours of data
        %   disp('creating the signal vector ...');
        %   tic
        %   X=sin([1:nsamp]/100)+sin([1:nsamp]*2)/10;
        %   X = [X' fliplr(X)']; % Two cahnnels
        %   toc
        %   filter=struct('order',2,'ripple',0.1,'stopBand',40,'low',600,'high',6000);
        %   [b,a]=ellip(filter.order,filter.ripple,filter.stopBand,[filter.low filter.high]./samplinghz);
        %   disp('performing filtfilt ...');
        %   tic
        %   Xff=filtfilt(b,a,X);
        %   toc
        %   disp('performing filtfiltpiecemeal ...');
        %   tic
        %   chunksamples=round(samplinghz*3600/4); % cut in 1.5 minute intervals
        %   marginsamples=max(samplinghz,round(samplinghz/filter.low*100));
        %   Xffpm=filtfiltpiecemeal(b,a,X,chunksamples,marginsamples);
        %   toc
        %   disp(['sum(abs(Xff-Xffpm))=' num2str(sum(abs(Xff-Xffpm)))]);
        %
        % See also: filtfilt
        %
        % Jacob, Jan-2012



        % input checks
        if rem(marginSamples,1)~=0
            error('marginsamples must be an integer');
        end
        if rem(chunkSamples,1)~=0
            error('chunksamples must be an integer');
        end
        if marginSamples>chunkSamples
            error('marginsamples must be less than chunksamples');
        end
        % main course
        begin=1; % keeps track of position in X, steps with chunksamples
        nchunks=ceil(nrSamples/chunkSamples);

        if nchunks ==1
            chunkSamples = nrSamples; % If you don't do this, the nrSamples beyond the first chunk are not used because the loop below never gets to cycle==chunks. BK - Aug. 2015
        end
        if nchunks==0
            signal=filtfilt(b,a,signal);
        else
            for cycle=1:nchunks
                %disp(['chunk: ' num2str(cycle)]);
                if cycle==1
                    arridx2change=1:min(begin+chunkSamples-1,nrSamples);
                    arridx2filter=1:min(begin+chunkSamples+marginSamples-1,nrSamples);
                    thisCyclePreMargin=zeros(0,nrSignals);
                    filtout2paste=1:min(chunkSamples,numel(arridx2filter));
                    nextCyclePreMarginIdx=min(begin+chunkSamples-marginSamples,nrSamples):min(begin+chunkSamples-1,nrSamples);
                    nextCyclePreMargin=signal(nextCyclePreMarginIdx,:);
                elseif cycle<nchunks
                    arridx2change=begin:min(begin+chunkSamples-1,nrSamples);
                    arridx2filter=begin:min(begin+chunkSamples+marginSamples-1,nrSamples);
                    thisCyclePreMargin=nextCyclePreMargin;
                    filtout2paste=( 1:min(chunkSamples,numel(arridx2filter)) )+size(thisCyclePreMargin,1);
                    nextCyclePreMarginIdx=min(begin+chunkSamples-marginSamples,nrSamples):min(begin+chunkSamples-1,nrSamples);
                    nextCyclePreMargin=signal(nextCyclePreMarginIdx,:);
                elseif cycle==nchunks
                    arridx2change=begin:nrSamples;
                    arridx2filter=begin:nrSamples;
                    thisCyclePreMargin=nextCyclePreMargin;
                    filtout2paste=( 1:numel(arridx2filter) )+size(thisCyclePreMargin,1);
                end

                if pv.parms.useGPU
                    % Copy chunk to GPU then process there and pull back
                    % filtfilt needs 5* the amount of memory needed for a
                    % single chunk. For a small GPU (4GB) the overhead is
                    % large and this does not lead to  a spedup.
                    filtout=filtfilt(b,a,gpuArray([ thisCyclePreMargin; signal(arridx2filter,:) ]));
                    filtout = gather(filtout);
                else
                    filtout=filtfilt(b,a,[ thisCyclePreMargin; signal(arridx2filter,:) ]);
                end

                signal(arridx2change,:)=filtout(filtout2paste,:);
                begin=begin+chunkSamples;
            end
        end
    end

if makeColumn
    % Revert to honor the input orientation
    signal=signal';
end


end
