%{
# Mapping samples to trials for preprocessed signals 
-> ns.Experiment
-> ephys.PrepParm
trial  : int
--- 
sample :  longblob   # Samples from this Preprocessed set that correspond to this trial
nstime   : longblob  # Time in seconds on the Neurostim experiment clock
trialtime   : longblob    # Time in seconds relatve to the first sample in the trial.
%}
% This determines how each sample in a session maps to a trial and a time in
% a Neurostim experiment. A sample is assigned to a trial if it occurs after
% the first monitor frame in the trial and before the first monitor frame of the next trial. 
% (In other words the ITI is included at the *end* of each trial).
% 
% Because different .Preprocessed data sets for the same ns.Experiment could have different
% numbers of samples, this is computed per Preprocessed set. 
%
% This table is used to extract activity per trial, aligned to first frame. 
% 

classdef PreprocessedTrialmap < dj.Part % Manual because it is automatically created by Preprocessed
 properties (SetAccess = protected)
        master = ephys.Preprocessed;  % Part  table for the plugin
    end
    methods  (Access = protected)
        function makeTuples(~,~)
             %Handled by the parent class ephys.Preprocessed
        end

    end
end