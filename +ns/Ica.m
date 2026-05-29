%{
# Stores Independent Component Analysis results for a row in ns.C
-> ns.C
-> ns.IcaParm
---
winverse : longblob  # Inverse ICA weights (mixing matrix)
sphere   : longblob  # Sphering matrix
weights  : longblob  # ICA weights (unmixing)
channels : longblob  # Channels used for ICA
%}
%
% Currently implemented with EEGLAB (s
%
% See also ns.IcaParm, ephys.eeglab.preprocess
%
% BK - May 2026
classdef Ica < dj.Computed
    properties (Dependent)
        keySource
    end

    methods
        function v = get.keySource(~)
            % Restricted to ns.C tuples with the ctag specified in IcaParm
            v = (ns.C * ns.CParm) * proj(ns.IcaParm,'itag','ctag');
        end
    end

    methods (Access=protected)
        function makeTuples(tbl,key)
            
            icaParms = fetch1(ns.IcaParm & key,'parms');
            % Get an EEG struct from ns.C using the ctag specified in IcaParm
            EEG = ephys.eeglab.dataset(key,data=key.ctag);
            % Use the parms from IcaParm to run ICA on the EEG struct, which will add the ICA fields to the struct. The parms should be a struct with a field 'eeglab' that contains a struct of parameters to pass to pop_runica. For example:
            parms.eeglab.ica =  icaParms; % Match the field name expected by preprocess.m
            EEG = ephys.eeglab.preprocess(EEG,parms);
            
            tpl = key;
            tpl.winverse = EEG.icawinv;
            tpl.sphere = EEG.icasphere;
            tpl.weights = EEG.icaweights;
            tpl.channels = EEG.icachansind;
            insert(tbl,tpl);
        end
    end
end
