%{
#  Preprocessing instructions for Spikes deconvolution
stag         :  varchar(255)     # A  unique name for these preprocessing instructions
---
deconv      : blob          # struct containing deconvolution parameters 
calibration = NULL : blob          # struct containing calibration parameters
%}
% Currently only for mlspike purposes
%  deconv.fluorescence can be used to indicate which ns.C ctag contains the
%  fluorescence data (defaults to "fluorescence")
%
%  Use a convenience function to setup the deconv struct with defaults
%  from the literature. 
%  
% mlparms = sbx.SpikesParm.defaults("DENEUX16","gcamp6s") 
% mlparms.deconv.algo.nspikemax =4;  % Allow 4 spikes per bin (15 Hz)
% autocalpamrs = struct('amin',0.035,...
%                'amax',0.2,...
%                 'taumin',0.25,...
%                 'taumax',2,...
%                 'maxamp',4);  % A maxamp of 4 seems necessary in our data.
% autocalparms.nrRoi = 20  % Maximum number of ROI to use for autocalibration of a session
%
% insertIfNew(sbx.SpikesParm,struct('stag','deneux16','deconv',mlparms,'calibration',autocalparms);
%
% 
% SEE ALSO sbx.mlspikecompute, Mlspikecalibration

classdef SpikesParm < dj.Lookup

    methods (Static)
        % A function to setup default parameters
        function deconv = defaults(publication,indicator,fineTuned)
            arguments
                publication (1,1) string
                indicator (1,1) string
                fineTuned (1,1) logical =false
            end

            deconv = tps_mlspikes('par');
            switch (publication)
                case "RUPPRECHT25"
                    % Return parameters for Gcamp6s based on Rupprecht et al. 2025.
                    switch upper(indicator)
                        case "GCAMP6S"
                            deconv.hill = 1.84;
                            deconv.a = 0.113;  % Amplitude
                            deconv.tau  =1.87; % Decay tau in s
                            deconv.ton = 0.07; % Rise tau (t_on) in s
                            deconv.pnonlin = 0.1;
                            deconv.drift.parameter = 0.1;
                        case "GCAMP8S"
                            deconv.hill = 2.2;
                            deconv.a = 0.576;  % Amplitude
                            deconv.tau  =0.267; % Decay tau in s
                            deconv.ton = 0.00472; % Rise tau (t_on) in s
                            deconv.pnonlin = 0.1;
                            deconv.drift.parameter = 0.1;
                            if fineTuned
                                % Use tau and amplitude as estimated by Rupprecht
                                % with a grid search on ground truth data
                                deconv.a    = 1;  % Amplitude
                                deconv.tau  = 0.7; % Decay tau in s
                            end
                        otherwise
                            error("No %s defaults for %s",indicator,publication)
                    end
                case "DENEUX16"
                    switch upper(indicator)
                        case "GCAMP6S"                            
                            deconv.hill = 1; % Using polynomial pnonlin instead
                            deconv.a = 0.07;  % Amplitude
                            deconv.tau  =1.3; % Decay tau in s
                            deconv.ton = 0.02; % Rise tau (t_on) in s
                            deconv.pnonlin =  [0.73, -0.05];  % Deneux values for Gcamp6s;
                            deconv.drift.parameter = 0.1;
                        otherwise
                            error("No %s defaults for %s",indicator,publication)
                    end
                otherwise 
                     error("Unknown publication %s",publication)
            end
        end
    end
end