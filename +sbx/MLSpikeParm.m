%{
#  Preprocessing instructions for MLSpike autocalibration and deconvolution
stag         :  varchar(255)     # A  unique name for these preprocessing instructions
---
parms       : longblob          # struct containing all parameters used by MLSpike
%}
%
% EXAMPLE
%
%mlParms =struct;
%mlParms.restrict = 'pcell>0.75 AND radius> 2.185';
%mlParms.deconv = sbx.SpikesParm.defaults("DENEUX2012") % Defaults defined in the class
%mlParms.deconv.dographsummary = false;
% Note that the "true" values for CA indicators often do not match the ones
% that lead to the best results ; autocalibration should help to find the best 
% parametsrs. 
% mlParms.deconv.algo.nspikemax =4;  % Allow 4 spikes per bin (15 Hz)
% mlParms.autocalibration = struct('amin',0.035,...
%                'amax',0.2,...
%                 'taumin',0.25,...
%                 'taumax',2,...
%                 'maxamp',4);  % A maxamp of 4 seems necessary in our data.
% insertIfNew(sbx.MLSpikeParm,struct('stag','mlspikeautocal','parms',mlParms);
%
% SEE ALSO sbx.Spikes
classdef MLSpikeParm < dj.Lookup

    methods (Static)
       end
end