function [exp_tpl, ch_tpl] = read_noisyChannels(cTbl, parms)
% read noisyChannel struct from ns.C info
% parms.include
% or
% parms.exclude
arguments

    cTbl (1,1) ns.C
    parms (1,1) struct

end

c_tpl = fetch(cTbl, '*');

% Check if noisyChannels exist
if isfield(c_tpl.info, "noisyChannels")

    noisies = c_tpl.info.noisyChannels;
    fld_names = string(fields(noisies));
    fld_names = fld_names(startsWith(fld_names,"badBy"));

    if ~isempty(parms)

        if isfield(parms, "include") && ~isempty(parms.include)

            assert(~isfield(parms, "exclude") || isempty(parms.exclude), ...
                "read_noisyChannels: parms should contain either 'include' or 'exclude' as a field, not both.");

            assert(all(ismember(parms.include, fld_names)), ...
                "read_noisyChannels: parms.include contains an unknown noisy channel category. Valid categories are\n %s",sprintf('\t%s\n',fld_names))


            fld_names = parms.include;

        elseif isfield(parms, "exclude") && ~isempty(parms.exclude)

            assert(~isfield(parms, "include") || isempty(parms.include), ...
                "read_noisyChannels: parms should contain either 'include' or 'exclude' as a field, not both.");
            assert(all(ismember(parms.exclude, fld_names)), ...
                "read_noisyChannels: parms.exclude contains an unknown noisy channel category. Valid categories are\n %s",sprintf('\t%s\n',fld_names));
            
            fld_names = fld_names(~ismember(fld_names, parms.exclude));

        end

    end

    bad_channels = arrayfun(@(fld_name) gen.make_row(noisies.(fld_name)), fld_names, 'UniformOutput', false);
    bad_channels = num2cell([bad_channels{:}]);
else
    warning("No record of noisy channels were found in ns.C.");
    bad_channels = [];
end

ch_tpl =struct('start',[],'stop',[],'trial',[],'channel',bad_channels);
exp_tpl = struct('start',[],'stop',[],'trial',[]); % Nothing applicable to all channels
end