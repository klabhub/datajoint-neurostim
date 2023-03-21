%{
# Parameters used to extract trial aligned activity from the Preprocessed data
act : varchar(25) #  A unique name for these parameter settings.
---
modality : varchar(20)  # 'spikes', 'fluorescence', 'neuropil'                
start : decimal(10,3)    # Start time in s relative to first frame in the trial 
stop  : decimal(10,3)    # Stop time in s relative to first frame in the trial 
step  : decimal(10,3)    # Step time in s. 
interpolation : varchar(20) # Interpolation method; see timetable/synchronize for allowed values.(e.g., 'linear','nearest','mean')
%}

classdef ActivityParms <dj.Lookup
    methods (Static)
        function example
            tpl = struct('act','spikes-nn', 'modality','spikes','start',-0.25,'step',0.050,'stop',2,'interpolation','nearest');
            insert(sbx.ActivityParms,tpl);
        end
    end
end