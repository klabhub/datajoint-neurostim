%{ 
# Calibrated MLSpike parameters for a session.
-> sbx.Preprocessed
-> sbx.SpikesParm
---
quality : float             # Correlation between reconstructed F and F
amplitude : float           # Calibrated F per spike
sigma : float               # Calibrated noise level    
tau : float                 # Calibrated decay parameter
failed : smallint           # Number of ROIs where calibration failed.
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
            % amplitude,sigma, tau. 
            parms = fetch(sbx.SpikesParm&key,'deconv','calibration');
            assert(all(isfield(parms.calibration,["amin" "amax" "taumin" "taumax" "maxamp" "nrRoi"])),'The %s SpikesParm does not have the required calibration parameters\n',parms.stag)

            if isfield(parms.deconv,'fluorescence')
                % Identifies a specific ns.C row
                fluorescence = parms.deconv.fluorescence;
            else
                % Default name for F
                fluorescence = "fluorescence";
            end
            % Query the CChannel table for fluorescence time series 
            allF = (ns.C & struct('ctag',fluorescence) )* (ns.CChannel &  proj(sbx.PreprocessedRoi & key,'roi->channel'));
            nrRoi = count(allF);
            nrCalibrationRoi = min(nrRoi,parms.calibration.nrRoi);
            calCntr = 0;
            calResults = repmat(struct('quality',[],'tau',[],'a',[],'sigma',[]),[nrCalibrationRoi 1]);
            info = sbx.readInfoFile(fetch(ns.Experiment & key,'LIMIT 1'));
            parms.deconv.dt = 1./(fetch1(sbx.Preprocessed & key,'framerate')/info.nrPlanes); % Match dt to framerate
            parms.calibration = rmfield(parms.calibration,"nrRoi");
            parms.calibration.dt = parms.deconv.dt;
            for thisF = fetch(allF,'*',sprintf('ORDER BY rand() LIMIT %d',nrCalibrationRoi))'
                    calCntr =calCntr+1;
                    calResults(calCntr) =    sbx.mlspikecompute(thisF,parms.deconv,parms.calibration);                    
            end

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