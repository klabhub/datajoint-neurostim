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
            
            
            if ~isfield(pv.parms, 'mesh')
               pv.parms.mesh.cfg  = struct;
            end
            assert(isstruct(pv.parms.mesh.cfg), 'Missing mesh.cfg in parms');
            if isempty(pv.parms.mesh.cfg) || ~isfield(pv.parms.mesh.cfg, 'downsample')
                pv.parms.mesh.cfg.downsample =1;
            end


            assert(isfield(pv.parms, 'sourcemodel'), 'Missing sourcemodel in parms');
            assert(isfield(pv.parms.sourcemodel, 'cfg'), 'Missing sourcemodel.cfg in parms');
            assert(isfield(pv.parms.sourcemodel.cfg, 'method'), 'Missing sourcemodel.cfg.method in parms');
            

            switch pv.parms.sourcemodel.cfg.method
                case 'basedonmri'
                    % Supplement defaults if missing.
                    if ~isfield(pv.parms.sourcemodel.cfg, 'resolution')                        
                        pv.parms.sourcemodel.cfg.resolution = 8; % in mm
                    end
                    if ~isfield(pv.parms.sourcemodel.cfg, 'tight')                        
                        pv.parms.sourcemodel.cfg.tight = 'no';
                    end
                    if ~isfield(pv.parms.sourcemodel.cfg, 'movetocentroids')                        
                        pv.parms.sourcemodel.cfg.movetocentroids = 'yes';
                    end
                case 'basedonresolution'
                    % Supplement defaults if missing.
                    if ~isfield(pv.parms.sourcemodel.cfg, 'resolution')                        
                        pv.parms.sourcemodel.cfg.resolution = 10; % in mm
                    end
                case 'basedonmni'
                    assert(isfield(pv.parms.sourcemodel, 'template'), 'Missing sourcemodel.template in parms');
                    assert(isfield(pv.parms.sourcemodel, 'atlas'), 'Missing sourcemodel,atlas in parms');
                    assert(isfield(pv.parms, 'elec'), 'Missing elec (template file) in parms');
         
                otherwise
                    error('Invalid sourcemodel.cfg.method in parms');
            end         
        end
    end
end