%{
# Epoching parameters to segment and preprocess trial data. 
etag : varchar (32) # unique tag
---
ctag            : varchar(32)       # C data that can be epoched with these parms
dimension       : varchar(32)       # Condition from the dimension table
window          : tinyblob          # Start and stop time of the epoch.
channels =NULL  : blob              # Channels to include. Defaults to all in the ctag
align           : blob              # struct defining the align event
prepparms       : blob              # struct with preprocessing parameters
artparms        : blob              # struct array with artifact removal parameters
plgparms        : blob              # struct array with plugin based trial selection parameters
%}

% MOz Feb, 2025
classdef EpochParm < dj.Lookup & dj.DJInstance
    methods
        function insert(self, tuples, varargin)
           % Overload the insert method to do argument validation and setup
           % defaults
            pv = namedargs2cell(tuples);
            tuples = self.validate(pv{:});
            % Call the superclass function after validation
            insert@dj.Lookup(self, tuples, varargin{:})
        end
    end

    methods (Static, Access = protected)
        function pv = validate(pv)
            arguments
                pv.etag (1,1) string
                pv.ctag  (1,1) string
                pv.dimension (1,1) string
                pv.window (1,2) 
                pv.channels (1,:) {mustBeNumeric} = []
                pv.prepparms (1,1) struct  = struct('enable',false);
                pv.artparms  (1,1) struct  = struct('enable',false);
                pv.plgparms (1,1) struct  = struct('enable',false);
                pv.align (1,1) struct =struct([]);
            end  
            % Check validity and set defaults
            pv.artparms = mustBeArtParm(pv.artparms);
            pv.plgparms = mustBePlgParm(pv.plgparms);
            pv.prepparms = mustBePrepParm(pv.prepparms);


             % validate 'dimension' and 'plugin' exist in ns.Dimension
            dimTbl = ns.Dimension & struct('dimension',pv.dimension);
            assert(count(dimTbl), ...
                'Dimension table does not contain dimension value of "%s"', pv.dimension);
            cTbl = ns.C & struct('ctag',pv.ctag);
            assert(count(cTbl), ...
                'C table does not contain ctag value of "%s"', pv.ctag);

            if isempty(pv.align)
                % Default to the startTime of the plugin that defined the
                % dimension. 
                %  Check that there is only one plugin for this dimension
                G = proj(dimTbl, 'dimension');                         
                multipleOptions = aggr(G, dimTbl, 'count(distinct plugin)->n') & 'n>1';
                assert(count(multipleOptions)==0,"The %s dimension links to multiple plugins. Cannot pick a default align.",pv.dimension);
                % Setup the align struct
                plg = fetch1(dimTbl,'plugin','LIMIT 1');
                pv.align = struct('plugin',plg{1},'event','startTime');
            end
                          
        end
    end

end

function prep= mustBePrepParm(prep)
 %TODO
end
function plg = mustBePlgParm(plg)
 %TODO
end

function art = mustBeArtParm(art)
% Validate artifact-detection parameters stored in ns.EpochParm.
% The validator is called by an arguments-block validator as a positional
% struct, so the input must be declared as a struct rather than as a
% name-value structure.
arguments
    art (1,1) struct
end

% See prep.artifactDetection for descriptions.
defaults = struct;
defaults.enable = false;
defaults.amplitude_threshold_peak = 150;
defaults.amplitude_channelfrac = 0.5;
defaults.variance_z_threshold = 5;
defaults.hf_cutoff_hz = 50;
defaults.hf_z_threshold = 5;
defaults.correlation_z_threshold = 3;
defaults.criterion_channels = [];
defaults.exclude = [];
defaults.epoch_no = [];
defaults.ica = [];

supplied = art;
art = defaults;
suppliedNames = fieldnames(supplied);
notAllowed = setdiff(string(suppliedNames),[string(fieldnames(defaults)); "flat_threshold_sd"]);
assert(isempty(notAllowed),"artparms do not allow %s fields",strjoin(notAllowed));
for iName = 1:numel(suppliedNames)
    art.(suppliedNames{iName}) = supplied.(suppliedNames{iName});
end

if isfield(art,'enable')
    validateattributes(art.enable,{'logical'},{'scalar'},mfilename,'enable');
end
numericScalarFields = ["amplitude_threshold_peak" "varianze_z_threshold" ...
    "hf_cutoff_hz" "hf_z_threshold" "correlation_z_threshold"];
for name = numericScalarFields
    if isfield(art,name)
        validateattributes(art.(name),{'double'},{'scalar','real','finite'}, ...
            mfilename,char(name));
    end
end
if isfield(art,'flag_if_channel_noisy') && ~isempty(art.flag_if_channel_noisy)
    validateattributes(art.flag_if_channel_noisy,{'double'}, ...
        {'row','real','finite'},mfilename,'flag_if_channel_noisy');
end
% The .ica struct is passed to ns.Ica/clean as pv.
if isfield(art,'ica') && ~isempty(art.ica)
    assert(isstruct(art.ica) && all(isfield(art.ica,["itag" "ltag"])), ...
        'The ICA artifact parameter must specify itag and ltag.');
    if isfield(art.ica,'find')
        assert(isempty(art.ica.find) || ...
            (isstruct(art.ica.find) && all(ismember(fieldnames(art.ica.find),["op" "name" "threshold"]))), ...
            'The ICA find operation must specify value or threshold and optionally the op.');
    end
end
end



