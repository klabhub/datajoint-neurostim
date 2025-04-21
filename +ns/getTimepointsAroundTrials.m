function varargout = getTimepointsAroundTrials(varargin)

% Subsamples signal to the segments +- buffer_t (time) around trials
% [isTimepointAroundTrials, subsampled_ timepoints] = getTimepointsAroundTrials(c_tbl, buffer_t)
% [isTimepointAroundTrials, subsampled_ timepoints] = getTimepointsAroundTrials(exp_tbl, timepoints, buffer_t)
% Outputs:
%   isTimepointAroundTrials
%   subsampled_ timepoints

if nargout == 2 && isa(varargin{1}, "ns.C")
    c_tbl = varargin{1};
    exp_tbl = ns.Experiment & c_tbl;
    buffer_t = varargin{2};
    assert(count(c_tbl)==1, "The function only accepts a single entry in C table.")
    t = sampleTime(c_tbl);
elseif nargout == 3 && isa(varargin{1}, "ns.Experiment")
    exp_tbl = varargin{1};
    t = varargin{2};
    buffer_t = varargin{3};
else
    error("If the first input is ns.C second and only other input must be buffer_t. If first input is ns.Experiment, second input must be timepoints, and third input must be buffer_t.");
end

assert(count(exp_tbl)==1, "The function only accepts a single entry in Experiment table.");
trl_t = get(exp_tbl, 'cic','prm','firstFrame','atTrialTime',inf,'what','clocktime');
n_trl = length(trl_t);

isInSubsamp = ones(size(t))==0;

for ii = 1:n_trl

    buffer_winN = trl_t(ii) + [-1, 1]*buffer_t;

    isInSubsamp = gen.ifwithin(t, buffer_winN) | isInSubsamp;

end

varargout = cell(1, nargout);
varargout{1} = isInSubsamp;
if nargout > 1, varargout{2} = t(isInSubsamp); end

end