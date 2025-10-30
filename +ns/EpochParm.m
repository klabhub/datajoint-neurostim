%{

# Epoching parameters to segment trial data. All times are relative to the plugin time

etag : varchar(32) # unique tag
paradigm : varchar(32)
dimension : varchar(32) # Condition from the dimension table
---
plugin : varchar(32) # The plugin to be timelocked to from DimensionCondition table
epoch_win : BLOB
pv : BLOB # structure array containing the custom parameters

%}

% MOz Feb, 2025
classdef EpochParm < dj.Lookup & dj.DJInstance

    methods

        function insert(self, tuples, varargin)

            % validate 'dimension' and 'plugin' exist in ns.Dimension
            dimTbl = ns.Dimension & sprintf('dimension ="%s"', tuples.dimension);
            assert(count(dimTbl), ...
                'Dimension table does not contain dimension value of "%s"', tuples.dimension);
            assert(ismember(tuples.plugin, dimTbl{"plugin"}), ...
                'Dimension table does not contain plugin value of "%s"', tuples.plugin);
            pv = namedargs2cell(tuples.pv);
            tuples.pv = self.validate_pv(pv{:});

            % Inherited function after validation
            insert@dj.Lookup(self, tuples, varargin{:})

        end
    end

    methods (Static, Access = protected)

        function pv = validate_pv(pv)

            arguments

                pv.resample (1,1) {mustBeNumeric} = 0
                pv.resample_opts = {}
                pv.detrend (1,1) {mustBeLogicalOrNumeric} = 0
                pv.baseline (1,1) {mustBeLogicalOrNumeric} = 0
                pv.baseline_win {validateBaselineWin} = []
                pv.rereference (1,1) {mustBeLogicalOrNumeric} = 0
                pv.ref_opts = {}
                pv.artifact_parm {validateArtParm} = struct(fun = '', args = {})

            end

            assert(~pv.baseline || ~isempty(pv.baseline_win), ...
                "Provide a baseline window (baseline_win) in the EpochParm entry or set 'baseline' to false.")

       end



    end

end

function mustBeLogicalOrNumeric(x)

assert(isnumeric(x) & (x == 0 || x == 1) || islogical(x), "Must be logical or 0 or 1.");

end

function validateBaselineWin(x)

assert(isempty(x) | numel(x) == 2, "Must be either empty or [t_start, t_end]");

end

function validateArtParm(art)

assert(isstruct(art) && all(ismember(fieldnames(art), ["fun", "args"])), ...
    "artifact_parm must be a struct with fields fun (artifact finder function) and args (input cell array)");
assert( ...
    all(all(cellfun(@(x)isstring(x) | ischar(x) | isa(x, 'function_handle'), {art.fun}))), ...
    "Field 'fun' must contain a string or handle associated with artifact finder function.");

assert(all(arrayfun(@(x) iscell(x) | isstruct(x), {art.args})), "Field args must contain cell array.");
end