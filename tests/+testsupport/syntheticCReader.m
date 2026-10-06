function [signal,time,channelInfo,recordingInfo] = syntheticCReader(key,parms)
%SYNTHETICCREADER Deterministic continuous data reader for integration tests.
arguments
    key (1,1) struct %#ok<INUSA>
    parms (1,1) struct
end

nSamples = parms.nSamples;
nChannels = numel(parms.channels);
time = [parms.startTime parms.stopTime nSamples];
sampleIndex = (0:nSamples-1)';
signal = zeros(nSamples,nChannels,'single');
for channel = 1:nChannels
    signal(:,channel) = single(sampleIndex + 1000*channel);
end

channelInfo = repmat(struct('nr',0,'name',''),1,nChannels);
for channel = 1:nChannels
    channelInfo(channel).nr = parms.channels(channel);
    channelInfo(channel).name = sprintf('synthetic-%d',parms.channels(channel));
end
recordingInfo = struct('source','syntheticCReader','nSamples',nSamples);
end
