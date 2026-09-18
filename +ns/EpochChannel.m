%{
# Preprocessed and epoched data per channel and trial/epoch
-> ns.Epoch
channel : int       # Channel number
trial : int         # Trial number 
---
onset : float             # Time of the align event relative to trial start
signal     : longblob     # C data for a single channel, single trial
%}
classdef EpochChannel < dj.Part & dj.DJInstance & ns.cache   
    properties (SetAccess = protected)
        master = ns.Epoch
    end

    properties (Dependent)
        channels                % Channels contributing to this Epoch table            
    end


    methods (Access = protected)
        function [src]= getCacheQuery(o)
             % Determine the complete query/relvar for epochs. The cache
             % class expects this to have the following columns:
             % time - the time of the samples, relative to the align event
             % align - the name of the event to which the epoch is aligned
             % signal - the actual data for the epoch
             % onset - the time of the align event relative to trial start
             % condition - the condition for each trial
             % paradgmm - the paradigm for each trial             
             src = proj(o,'signal','onset') * proj(ns.Epoch,'time') * proj(ns.EpochParm,'align') * proj(ns.Experiment,'paradigm') * proj(ns.DimensionTrial,'name->condition');             
        end
    end

    methods      
        function insert(self,tuples,varargin)
            insert@ns.cache(self,tuples);
            insert@dj.Part(self,tuples,varargin{:});
        end

        function ch = get.channels(self)
            ch = self.unique('channel');
        end               
    end

   
   

    
   

end
