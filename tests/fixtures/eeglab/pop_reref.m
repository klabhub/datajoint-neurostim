function EEG = pop_reref(EEG, ~, varargin)
% Test-only shim for the EEGLAB pop_reref function.
if ~isfield(EEG, 'etc') || ~isfield(EEG.etc, 'rerefCalls')
    EEG.etc.rerefCalls = [];
end
EEG.etc.rerefCalls(end+1) = numel(EEG.etc.rerefCalls) + 1;
assert(isempty(varargin), 'The default rereference should not pass parameters.');
end
