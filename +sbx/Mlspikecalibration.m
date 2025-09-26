%{ 
# Calibrated MLSpike parameters for a session.
-> sbx.Preprocessed
-> sbx.SpikesParm
---
quality = NULL: float             # Correlation between reconstructed F and F
a       = NULL: float           # Calibrated F per spike
sigma   = NULL: float               # Calibrated noise level    
tau     = NULL: float                 # Calibrated decay parameter
failed  : smallint           # Number of ROIs where calibration failed.
%}
classdef Mlspikecalibration < dj.Computed
    properties (Dependent)
        keySource
    end

    methods 
        function v = get.keySource(tbl)
            v = sbx.Preprocessed * (sbx.SpikesParm & 'calibration IS NOT NULL ');
        end
    end
    methods (Access = protected)
        function makeTuples(tbl,key)
            % Select a subset of the roi in the session, run
            % autocalibration and then pick the median of the calibrated 
            % amplitude,sigma, tau and store this in the table for later use 
            % when populating the sbx.Spikes table.

            calResults = sbx.mlspike(key,calibration =true);
            key.tau = mean([calResults.tau],"omitmissing");
            key.sigma  =mean([calResults.sigma],"omitmissing");
            key.a = mean([calResults.a],"omitmissing");
            key.quality = mean([calResults.quality],"omitmissing");
            key.failed = sum(isnan([calResults.quality]));
            % Store the calibration results             
            insert(tbl, key);

        end
    end
end