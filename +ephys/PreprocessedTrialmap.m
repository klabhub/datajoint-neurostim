%{
# Mapping samples to trials for preprocessed signals 
-> ephys.Preprocessed
trial  : int
--- 
startsample :  int # First sample that belongs to this trial
stopsample  : int  # Last sample that belongs to this trial
trialstart: float # Time in seconds of the first sample (on the Neurostim experiment clock)
sampleduration  : float # Duration of a single sample
%}
% This contains the information how each sample in a session (stored in 
% the PreprocessedChannel table) maps to a trial and time in the
% experiment.
% A sample is assigned to a trial if it occurs after
% the first monitor frame in the trial and before the first monitor frame of the next trial. 
% (In other words the ITI is included at the *end* of each trial).
% 
% Because different .Preprocessed data sets for the same ns.Experiment could have different
% numbers of samples, this is computed per Preprocessed set. 
%
% This table is used to extract activity per trial, aligned to first frame. 
% See also ephys.Preprocessed.get
classdef PreprocessedTrialmap < dj.Part 
 properties (SetAccess = protected)
        master = ephys.Preprocessed;  
    end
    methods  (Access = protected)
        function makeTuples(~,~)
             %Handled by the parent class ephys.Preprocessed
        end

    end
end