%{
# Transformed epoched data per channel 
-> ns.Tepoch
channel : int       # Channel number
trial : int         # Trial number 
---
y : longblob         # Data 
%}
classdef TepochChannel < dj.Part & dj.DJInstance & ns.cache    
    properties (SetAccess = protected)
        master = ns.Tepoch
    end 

    properties (Dependent)
        channels                % Channels contributing to this Epoch table
        samplingRate 
    end

    methods (Access = protected)
        function src = getCacheQuery(o)
            src = o*ns.Tepoch* proj(getCacheQuery(ns.EpochChannel & o),'condition','paradigm','align');
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

   