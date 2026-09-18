%{
# Transformed epoched data per channel and per trial.
-> ns.Tepoch
channel : int       # Channel number - can be zero to represent an average of channels
trial : int         # Trial number - can represent the first trial in a set of trials with the same condition
group = "_" : varchar(32) # A name that identifies a group of trials or channels that were averaged together. If not averaged, this is '_'
---
signal : longblob         # (Transformed) Data 
nrtrials = 1 : int       # Number of trials (if averaged)
nrchannels = 1 : int      # Number of channels (if averaged)
%}
classdef TepochChannel < dj.Part & dj.DJInstance & ns.cache    
    properties (SetAccess = protected)
        master = ns.Tepoch
    end 

    properties (Dependent)
        channels                % Channels contributing to this TEpoch table        
    end

    methods (Access = protected)
        function src = getCacheQuery(o)
        % Determine the complete query/relvar for tepochs. The cache
        % class expects this to have the following columns:
        % time - the time of the samples, relative to the align event
        % align - the name of the event to which the epoch is aligned
        % signal - the actual data for the epoch
        % onset - the time of the align event relative to trial start
        % If the Tepoch is grouped/averaged, condition column represents the name of the group of averaged trials and/or channels.
        % If the Tepoch is not grouped/averaged, condition column represents the condition for each trial (same as the condition of the related epoch)        
            groups = string(fetchn(o,'group'));
            if all(groups=="_") 
                % No groups, use condition from dimension
                src = proj(o,'signal') * proj(ns.Tepoch,'x','independent','dependent') * proj(ns.Epoch,'time') * proj(ns.EpochParm,'align') * proj(ns.Experiment,'paradigm') * proj(ns.DimensionTrial,'name->condition');                             
            else
                % Use groups as conditions
                src = proj(o,'group->condition','dependent','signal') * proj(ns.Tepoch *ns.Epoch * ns.EpochParm * proj(ns.Experiment,'paradigm'),'independent','x','time','align','paradigm');
            end
        end
    end

    methods      
        function ch = get.channels(self)
            ch = fetch(self, 'channel');
            ch = unique([ch(:).channel]');
        end
        
    
    end
end

   