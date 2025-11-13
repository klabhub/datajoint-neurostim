function [isTimepointAroundTrials, subsampled_timepoints] = getTimepointsAroundTrials(nsTbl, buffer, options)
% Subsamples signal to the segments defined by buffer around trials.
% [isAround, t_sub] = getTimepointsAroundTrials(exp_tbl, buffer, timepoints=t, plugin=plugin1)
%
% --- Inputs ---
%   nsTbl       : (1,1) ns.C OR ns.Experiment
%   buffer      : (1,:) double. 
%                 If scalar: [start-buffer, stop+buffer]
%                 If 2-element: [start-buffer(1), stop+buffer(2)]
%   timepoints  : (1,:) double = missing. Required if nsTbl is ns.Experiment.
%   plugin      : text = 'cic'. The plugin name used in the get() query.

arguments
    nsTbl {mustBeA(nsTbl, ["ns.C", "ns.Experiment"])}
    buffer (1,:) double {mustBeNonempty}
    options.timepoints (1,:) double = missing
    options.plugin {mustBeText} = 'cic'
end

% --- 1. Process Buffer ---
if isscalar(buffer)
    % Symmetric buffer
    buf_pre = - buffer;
    buf_post = buffer;
elseif numel(buffer) == 2
    % Asymmetric buffer: buffer(1) before start, buffer(2) after stop
    buf_pre = buffer(1);
    buf_post = buffer(2);
else
    error("Buffer must be a scalar or a 2-element vector.");
end

% --- 2. Process Table Type & Timepoints ---
if isa(nsTbl, "ns.C")
    % Syntax 1: ns.C
    assert(count(nsTbl)==1, "The function only accepts a single entry in C table.");
    
    % Derive Experiment and Time from C table
    exp_tbl = ns.Experiment & nsTbl;
    t = sampleTime(nsTbl);
    
    % Ensure t is (1,:) as requested by the timepoints signature
    if size(t,1) > 1 && size(t,2) == 1
        t = t'; 
    end

elseif isa(nsTbl, "ns.Experiment")
    % Syntax 2: ns.Experiment
    if ismissing(options.timepoints)
        error("If the first input is ns.Experiment, 'timepoints' (arg3) must be provided.");
    end
    
    exp_tbl = nsTbl;
    t = options.timepoints;
end

assert(count(exp_tbl)==1, "The function only accepts a single entry in Experiment table.");

% --- 3. Get Trial Timings ---
% Using the 'plugin' variable in the get query
trl_start_t = get(exp_tbl, options.plugin, 'prm', 'firstFrame', 'atTrialTime', inf, 'what', 'clocktime');
trl_stop_t  = get(exp_tbl, options.plugin, 'prm', 'trialStopTime', 'atTrialTime', inf, 'what', 'clocktime');

n_trl = length(trl_start_t);
isInSubsamp = false(size(t)); 

% --- 4. Loop & Threshold ---
for ii = 1:n_trl
    % Apply buffer: Start - buf_pre, Stop + buf_post
    curr_win = [trl_start_t(ii) + buf_pre, trl_stop_t(ii) + buf_post];
    
    % Accumulate logical mask
    isInSubsamp = do.ifwithin(t, curr_win) | isInSubsamp;
end

% --- 5. Outputs ---
isTimepointAroundTrials = isInSubsamp;
subsampled_timepoints = t(isInSubsamp);

end