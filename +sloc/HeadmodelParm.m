%{
# Parameters used to construct headmodels for source localization with fieldtrip.
hmtag : varchar (32) # unique tag
---
description: varchar(512) # A brief description of this model
parms : longblob # The parameters for the headmodel
%}

classdef  HeadmodelParm < dj.Lookup & dj.DJInstance
    methods
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
                pv.hmtag (1,1) string
                pv.description (1,1) string
                pv.parms  (1,1) struct 
            end
            
            defaults.mesh.cfg.downsample = 2; % This .cfg will be passed to ft_prepare_mesh
            defaults.sourcemodel.cfg.resolution = 8;  % This .cfg will be passed to ft_prepare_sourcemodel
            defaults.sourcemodel.cfg.method = 'basedonmni';
            defaults.sourcemodel.template = 'standard_sourcemodel3d5mm'; % Only needed/used for basedonmni
            defaults.sourcemodel.atlas= 'brainnetome/BNA_MPM_thr25_1.25mm.nii';% Only needed/used for basedonmni

            pv.parms= mergeDefaults(pv.parms,defaults);
           
        end
    end
end