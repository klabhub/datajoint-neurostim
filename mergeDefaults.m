function parms = mergeDefaults(parms, defaults)
% Recursively fill missing fields in PARMS from DEFAULTS.
arguments
    parms (1,1) struct
    defaults (1,1) struct
end

defaultFields = fieldnames(defaults);
for iField = 1:numel(defaultFields)
    name = defaultFields{iField};
    if ~isfield(parms, name)
        parms.(name) = defaults.(name);
    elseif isstruct(parms.(name)) && ...
            isscalar(parms.(name)) && isstruct(defaults.(name)) && ...
            isscalar(defaults.(name))
        parms.(name) = sloc.SourceParm.mergeDefaults( ...
            parms.(name), defaults.(name));
    end
end
end