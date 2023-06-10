%{
# Map samples in a preprocessed data set to trials in Experiments and time.
->ephys.Preprocessed   # Which preprocessed data set does this apply to
trial  : int
--- 
sample :  blob   # Samples from this Preprocessed set that correspond to this trial
nstime   : blob  # Time in seconds on the Neurostim experiment clock
trialtime   : blob    # Time in seconds relatve to the first sample in the trial.
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
% BK - June 2023
classdef PreprocessedTrialmap < dj.Part
    properties (SetAccess = protected)
        master = ephys.Preprocessed
    end 
    
end