%{
# Preprocessed signals per electrode and experiment 
-> ns.Experiment
-> ephys.PrepParm
channel : int # The channel that recorded the signal
---
signal : longblob # The preprocessed signal for each sample in the experiment
info = NULL : longblob # Any additional information on this preprocessed channel
%}
%
% This is a generic table that can store preprocessed data for a range  of
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
classdef Preprocessed < dj.Imported
    properties (Dependent)
        %    keySource
    end
    methods


    end

    methods (Access=protected)
        function makeTuples(tbl,key)
            % Evaluate the prep function for the specified parms

            parms = fetch(ephys.PrepParm & key,'*');
            % The (user-provided) prep funciton has to return the signal
            % as [nrSamples, nrChannels]  and the channels [nrChannels 1]
            [signal,channel,time,info] = feval(parms.fun,key,parms.parms);
channel =1;
            tpl = mergestruct(key,...
                struct('signal',num2cell(signal,1),...
                'channel',num2cell(channel(:))));
            if ~isempty(info)
                if iscell(info) &&  numel(info)==numel(channel)
                    % Add it
                    [tpl.channel] = deal(info{:});
                else
                    error('The channel info returned by the preprocessing function %s does not match the number of channels.',parms.fun)
                end
            end

            insert(tbl,tpl)

            %% Map samples to trials
            
            nrTrials = fetch1(ns.Experiment &key,'trials');
            % Events in neurostim are aligned to firstFrame; we do the same
            % for the analog data
            prms =  get(ns.Experiment & key,'cic');
            trialStartTime = (prms.cic.firstFrameNsTime/1000);
            % Create a cell array of logical addresses
            stay = cell(1,nrTrials);
            for tr=1:nrTrials
                stay{tr} = time >=trialStartTime(tr);
                if tr<nrTrials
                    stay{tr} = stay{tr} & time < trialStartTime(tr+1);
                end
            end
            thisSamples = 1:numel(time);
            % Samples per trial
            samples = cellfun(@(x)(thisSamples(x)),stay,'uni',false)';
            % Neurostim clock time per trial
            nsTimes = cellfun(@(x)(time(x)),stay,'uni',false)';
            % Tiem in the trial per trial (aligned to firstFrame in
            % neurostim)
            trialTimes= cellfun(@(x,y) (time(x)-y),stay,num2cell(trialStartTime)','uni',false)';
            % Create tuples and insert.
            tpl = struct('subject',key.subject,...
                'session_date',key.session_date,...
                'prep',key.prep,...
                'starttime',key.starttime,...
                'trial',num2cell(1:nrTrials)',...
                'sample',samples, ...
                'nstime',nsTimes,...
                'trialtime',trialTimes);
            insert(ephys.PreprocessedTrialmap,tpl);
        end
    end

end