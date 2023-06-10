%{
# Preprocessed signals per electrode and experiment 
-> ns.Experiment
-> ephys.PrepParm
--- 
%}
classdef Preprocessed < dj.Imported
    properties (Dependent)
        %    keySource
    end
    methods


    end

    methods (Access=protected)
        function makeTuples(tbl,key)
        
            insert(tbl,key)
            % Evaluate the prep function for the specified parms
            parms = fetch(ephys.PrepParm & key,'*');
            % The (user-provided) prep funciton has to return the signal
            % as [nrSamples, nrChannels]  and the channels [nrChannels 1]
            [signal,channel,time,info] = feval(parms.fun,key,parms.parms);
            [nrSamples,nrChannels] = size(signal);
            channel =1;
            assert(nrSamples==numel(time),'The number of rows in the preprocessed signal does not match the number of time points ')
            assert(nrChannels==numel(channel),'The number of columns in the preprocessed signal does not match the number of channels')
            assert(isempty(info) || (iscell(info )&& numel(info)==nrChannels),'The info rerturned by preprocessing deos not match the number of channels')
            
            channelsTpl = mergestruct(key,...
                       struct('signal',num2cell(signal,1),...
                              'channel',num2cell(channel(:))));
            if ~isempty(info)
                % Add it
                [channelsTpl.channel] = deal(info{:});
            end

            insert(ephys.PreprocessedChannel,channelsTpl);

           
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
            thisSamples = 1:nrSamples;
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