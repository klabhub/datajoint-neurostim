%{
# Parameters for SSVEP source localization.
stag : varchar(10) # Unique tag
---
description: varchar(128)  # Description 
parms : longblob # A struct with parameters.
%}

classdef SourceParm  < dj.Lookup
    
    methods (Access = public)
         function insert(self, tuples, varargin)
           % Overload the insert method to do argument validation and setup
           % defaults before insertion
            pv = namedargs2cell(tuples);
            tuples = self.validate(pv{:});
            % Call the superclass function after validation
            insert@dj.Lookup(self,makeMymSafe( tuples), varargin{:})
        end

    end 

     methods (Static, Access = protected)
        function pv = validate(pv)
            arguments
                pv.stag (1,1) string
                pv.description (1,1) string = ""
                pv.parms  (1,1) struct 
            end
            if upper(pv.parms.mode) =="BASELINE"
                % Need baseWindow  to do baseline correction of the sources
                requiredFields = ["baseWindow" "responseWindow" "mode","epoch"];
            elseif upper(pv.parms.mode) == "SNR" 
                requiredFields = ["responseWindow" "mode","epoch"];
            else
                error('Unknown mode %s',pv.parms.mode)
            end
            assert(all(isfield(pv.parms,requiredFields)),"Parms needs %s",strjoin(requiredFields));
            switch (pv.parms.mode)
                case {"BASELINE", "SNR"}
                    defaults.harmonic = 1;
                    defaults.signalHalfWidth = 2;
                    defaults.noiseHalfWidth  = 4;
                    defaults.freq.method      = 'mtmfft';                    
                    defaults.freq.taper       = 'rectwin';
                    defaults.freq.keeptrials  = 'yes';
                    defaults.freq.tapsmofrq   = 2;
                    defaults.freq.pad= 'nextpow2';

                    defaults.source.method                  = 'dics';
                    defaults.source.dics.projectnoise       = 'yes';
                    defaults.source.dics.lambda             = '5%';
                    defaults.source.dics.keepfilter         = 'yes';
                    defaults.source.dics.realfilter         = 'yes';
                    defaults.source.dics.fixedori           = 'yes';

                    pv.parms = mergeDefaults(pv.parms, defaults);

                    assert(isscalar(pv.parms.signalHalfWidth) && ...
                        isfinite(pv.parms.signalHalfWidth) && ...
                        pv.parms.signalHalfWidth > 0, ...
                        "signalHalfWidth must be a positive finite scalar.");
                    assert(isscalar(pv.parms.noiseHalfWidth) && ...
                        isfinite(pv.parms.noiseHalfWidth) && ...
                        pv.parms.noiseHalfWidth > pv.parms.signalHalfWidth, ...
                        "noiseHalfWidth must be greater than signalHalfWidth.");

                otherwise
                    error('Unknown source localization mode %s.', pv.parms.mode);
            end

        end

        
    end
end
