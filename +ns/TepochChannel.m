%{
# Transformed epoched data per channel and per trial.
-> ns.Tepoch
channel : int       # Channel number - can be zero to represent an average of channels
trial : int         # Trial number - can represent the first trial in a set of trials with the same condition
---
signal : longblob         # (Transformed) Data 
nrtrials = 1 : int       # Number of trials (if averaged)
nrchannels = 1 : int      # Number of channels (if averaged)
%
% NOTE:
% If the tepoch is created as an average over channels (TepocParm.average= 'channel') then
% the channel number is set to zero and nrchannels reflects the 
% number of distinct channels in the set that was averaged.
% 
% if the tepoch is created as an average over trials (tepochParm.average = 'trial')
% then the trial number represents the first trial for the condition and
% nrtrials is the number of trials in the average. This allows a join like:
% ns.Tepoch*proj(ns.DimensionTrial,'name->condition') to get the condition
% for each Tepoch.
%
%}
classdef TepochChannel < dj.Part & dj.DJInstance & ns.cache    
    properties (SetAccess = protected)
        master = ns.Tepoch
    end 

    properties (Dependent)
        channels                % Channels contributing to this TEpoch table
        samplingRate 
    end

    methods (Access = protected)
        function src = getCacheQuery(o)
            src = o*ns.Tepoch* proj(getCacheQuery(ns.EpochChannel & proj(o)),'paradigm','align');
        end
    end

    methods      
        function ch = get.channels(self)
            ch = fetch(self, 'channel');
            ch = unique([ch(:).channel]');
        end
        
        function v = get.samplingRate(self)
            e = ns.Epoch &self;
            v = uniquetol(e.samplingRate,0.01);
        end   
    end
end

   