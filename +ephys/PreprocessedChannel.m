%{
# Preprocessed signals per channel 
-> ephys.Preprocessed
channel : int # The channel that recorded the signal
---
signal : longblob # The preprocessed signal for each sample in the experiment
info = NULL : longblob # Any additional information on this preprocessed channel
%}
%
% This is a generic table that stores preprocessed data for a range  of
% electrophysiological recordings. For a specific signal the user provides
% the name of a function that reads data, preprocesses them and the returns
% values that are stored in the preprocessed table. Specifically, this
% preprocessing funciton has the following prototype:
% [signal,channel,time,info] = myprep(key, parms)
% INPUT
% key   -  The tuple that defines the source data (i.e. a row in the join of ns.Experiment
%               and ephys.PrepParm). The user function uses this key
%               to laod the appropriate file  (e.g. from ns.File & key),
%
% parms  - A struct with parameters (from (ephys.PrepParm &key)). The user
%           function uses this to determine what kind of preprocessing to do.
% OUTPUT
% signal -  A matrix with samples along the rows and channels along
%           columns.  These are all samples obtained during the experiment.
% channel  - A row vector with the channel number corresponding to each
%               column in signal.
% time -  The time at which each sample was obtained. This is a [nrSamples 1]
%               column vector. Note that the clock for these times must
%               match the neurostim clock. The user function likely
%               includes some code to matchup the clock of the data
%               acquisition device with some event in Neurostim.
% info - This is an optional cell array of size [1 nrChannels] each cell
%           provides some additional information on the channel that is stored in the
%           info field of the Preprocessed table. For example, this coudl contain the
%           hardware filtering parameters of the channel as read from the raw data
%           file.
%
% EXAMPLE:
% The ephys.ripple.prep functon shows a complete implementation of a
% preprocessing function that handles MUAE, LFP, and EEG recordings with
% the Ripple Grapevine system.
%
% BK - June 2023

classdef PreprocessedChannel < dj.Part
     properties (SetAccess = protected)
        master = ephys.Preprocessed;  % Part  table for the plugin
    end

    methods  (Access = protected)
        function makeTuples(~,~)
             %Handled by the parent class ephys.Preprocessed
        end

    end
end