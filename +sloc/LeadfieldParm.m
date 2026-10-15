%{
# Parameters used to construct leadfields for source localization with fieldtrip.
lftag : varchar (32) # unique tag
---
parms: longblob # The parameters for the leadfield computation
description: varchar(512) # A brief description of this leadfield computation
%}
classdef  LeadfieldParm < dj.Lookup & dj.DJInstance    
    methods
        function insert(self, tuples, varargin)
           % Overload the insert method to do argument validation and setup
           % defaults before insertion
            pv = namedargs2cell(tuples);
            tuples =validate(pv{:});
            % Call the superclass function after validation
            insert@dj.Lookup(self,makeMymSafe( tuples), varargin{:})
        end
    end
end

function pv = validate(pv)
arguments
    pv.lftag (1,1) string
    pv.description (1,1) string = "No description provided"
    pv.parms  (1,1) struct
end
parms = namedargs2cell(pv.parms);
pv.parms = mustBeLeadfieldParm(parms{:});

end

function  pv = mustBeLeadfieldParm(pv)
    %.mode is used to determine whether to loop over electrodes (e2d) or dipoles (d2e) or determine which one is faster (auto)
    % If a parallel pool is available, these modes will use it (see nsParPool) 
    % To force using the FieldTrip method, set mode to 'ft' (no parallelization, always per electrode, slow)
    %.cfg is passed to ft_prepare_leadfield and should only contain valid parameters for that function.
    % See ft_prepare_leadfield for details on the parameters.
    % Set explicit defaults if not specified.
    arguments
        pv.mode (1,1) string {mustBeMember(pv.mode,["auto", "e2d","d2e" "ft"])} = "auto"; 
        pv.cfg (1,1) struct  = struct('normalize','yes', 'normalizeparam', 0.5, 'weight', 1, 'reducerank', 'no', 'backproject', 'yes');
    end 

    % Set explicit defaults for the cfg fields if they are not specified.
    if ~isfield(pv.cfg,'normalize')
        pv.cfg.normalize = 'yes';
    end
    if ~isfield(pv.cfg,'normalizeparam')
        pv.cfg.normalizeparam = 0.5;
    end
    if ~isfield(pv.cfg,'weight')
        pv.cfg.weight = 1;
    end
    if ~isfield(pv.cfg,'reducerank')
        pv.cfg.reducerank = 'no';
    end
    if ~isfield(pv.cfg,'backproject')
        pv.cfg.backproject = 'yes';
    end

    allowedFields = [
        "normalize"
        "normalizeparam"
        "weight"
        "reducerank"
        "backproject"
        ];

    if any(~ismember(fieldnames(pv.cfg), allowedFields))
        error('Invalid field in parms. Allowed fields are: %s', strjoin(allowedFields, ', '));
    end
    
    
end