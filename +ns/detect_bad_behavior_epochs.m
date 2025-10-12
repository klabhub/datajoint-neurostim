function varargout = detect_bad_behavior_epochs(varargin)

warning(['Deprecation Warning!!\n\t ns.detect_bad_behavior_epochs is now' ...
    ' ns.prep.retrieve_bad_behavior_flags and will be removed in the future.']);

varargout = cell(1, nargout);
[varargout{:}] = ns.prep.retrieve_bad_behavior_flags(varargin{:});

end