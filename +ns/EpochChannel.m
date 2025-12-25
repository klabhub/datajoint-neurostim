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