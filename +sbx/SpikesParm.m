%{
#  Preprocessing instructions for Spikes deconvolution
stag         :  varchar(255)     # A  unique name for these preprocessing instructions
---
parms       : longblob          # struct containing all parameters used by packages such as mlspike.
%}
%
% EXAMPLE
%
%mlParms =struct;
%mlParms.restrict = 'pcell>0.75 AND radius> 2.185';
%mlParms.deconv = tps_mlspikes('par'); % Defaults
%mlParms.deconv.dographsummary = false;
%mlParms.deconv.pnonlin = [0.85 -0.006];  % Deneux values for Gcamp6s
%mlParms.deconv.drift.parameter = 0.01; % 
% a , tau, and sigma are determined by autocalibration, but if that fails
% (no isolated events) then these defaults weill be used:
%mlParms.deconv.a = 0.07;
% mlParms.deconv.tau  = 0.4; % Gcamp6s 1.3s but 0.8 is a typical value in
%                               our recordings
% mlParms.deconv.algo.nspikemax =4;  % Allow 4 spikes per bin (15 Hz)
% mlParms.autocalibration = struct('amin',0.035,...
%                'amax',0.2,...
%                 'taumin',0.25,...
%                 'taumax',2,...
%                 'maxamp',4);  % A maxamp of 4 seems necessary in our data.
% insertIfNew(sbx.SpikesParm,struct('stag','mlspikeautocal','parms',mlParms);
%
% SEE ALSO sbx.Spikes
classdef SpikesParm < dj.Lookup

end