%{
# Preprocessed and epoched data per channel and trial/epoch
-> ns.Epoch
channel : int       # Channel number
trial : int         # Trial number 
---
condition : varchar(128)  # Name of the condition to which this trial belongs
onset : float             # Time of the align event relative to trial start
signal : longblob         # C data for a single channel, single trial
%}
classdef EpochChannel < dj.Part & dj.DJInstance & ns.cache   
    properties (SetAccess = protected)
        master = ns.Epoch
    end

    properties (Dependent)
        channels                % Channels contributing to this Epoch table    
        samplingRate
    end

    methods (Access=public)
       

    end
    methods (Access = protected)
        function [src]= getCacheQuery(o)
            
             src = o * ns.Epoch * ns.EpochParm * proj(ns.Experiment,'paradigm');
             % Safety check; time and align should match for all rows
             % in the table. 
             preFetch = fetchtable(src,'time','align');
             assert(isscalar(unique(preFetch.time(:,3))),'Rows of the EpochChannel table must have the same numbers of samples.');
             assert(isscalar(unique({preFetch.align.plugin})),'Rows of the EpochChannel should be aligned to the same plugin.');
             assert(isscalar(unique({preFetch.align.event})),'Rows of the EpochChannel should be aligned to the same event.');          
             src = proj(src,'time','align','signal','condition','paradigm');
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