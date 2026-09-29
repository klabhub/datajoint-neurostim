function art = mustBeArtParm(art)
% Validate artifact-detection parameters stored in ns.EpochParm.
% The validator is called by an arguments-block validator as a positional
% struct, so the input must be declared as a struct rather than as a
% name-value structure.
    arguments
        art (1,1) struct
    end

    % See prep.artifactDetection for descriptiosn.
    defaults = struct;
    defaults.enable = false;
    defaults.amplitude_threshold_peak = 1;
    defaults.amplitude_channelfrac = 1;
    defaults.variance_z_threshold = 5;
    defaults.hf_cutoff_hz = 50;
    defaults.hf_z_threshold = 5;
    defaults.correlation_z_threshold = 3;
    defaults.criterion_channels = [];
    defaults.exclude = [];
    defaults.epoch_no = [];

    defaults.ica = struct(name = string.empty,threshold=NaN,op=function_handle.empty);

    supplied = art;
    art = defaults;
    suppliedNames = fieldnames(supplied);
    notAllowed = setdiff(string(suppliedNames),string(fieldnames(defaults)));
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