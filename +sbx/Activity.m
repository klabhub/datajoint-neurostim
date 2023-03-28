%{
# The activity of an ROI in one trial of an experiment
-> ns.Experiment
-> sbx.Roi
-> sbx.ActivityParms
trial :smallint
---
activity : longblob   # Activity trace per trial
%}
% This table stores activity traces per Experiment, per ROI and per trial.
% By using a specific ActitivityParms row, this activity can be spikes, or
% fluorescence, it can span different time periods in a trial, or sampled differently.
%
% The actual activity data (e.g., spikes) duplicate what is already stored
% in the sbx.ROI table and can also be retrieved from there using sbx.ROI.get.
%
% This table is populated by a specific call to sbx.ROI.get (as determined by
% the ActiviityParms; this saves time for subsequent analyses.
%
classdef Activity < dj.Computed

    methods
        function m = tuning(tbl,expts, conditions, neurons, pv)
            arguments
                tbl (1,1) sbx.Activity  
                expts
                conditions 
                neurons 
                pv.polar (1,1) logical = false
                pv.act {mustBeText} = 'spikes';
                pv.start = seconds(1)
                pv.stop = seconds(3)
            end
            if iscellstr(conditions) 
                x = get(ns.Experiment & expts,conditions{1},'prm',conditions{2},'attrialtime',0);
                [uX,~,xIx] = unique(x);
                [nrTrials,nrExpts] = size(x);
                nrX = numel(uX);
            elseif isa(conditions,'ns.Condition')
            
            end
            
            % Get the spikes Activity data for these roi
            [spk,time] = get(tbl & expts & neurons, pv.act);
           [nrTimePoints,nrTrials, nrNeurons, nrExpts] =size(spk);

            window = isbetween(time,pv.start,pv.stop);
            meanResponseInWindow = squeeze(mean(spk(window,:,:,:),1,"omitnan")); % Average over window
            meanResponseInWindow  = permute(meanResponseInWindow ,[1 3 2]); % [nrTrials nrExpts nrRois]
            meanResponseInWindow  = reshape(meanResponseInWindow ,nrTrials,[]); % [nrTrials nrRois]
            m = nan(nrX,nrNeurons);
            se = nan(nrX,nrNeurons);
            for c=1:nrX
                stay = x==uX(c);
                trialsPerCondition = sum(stay);
                m(c,:) = mean(meanResponseInWindow(stay,:),1,"omitnan");
                se(c,:) = std(meanResponseInWindow(stay,:),0,1,"omitnan")/sqrt(trialsPerCondition);
            end
           m = m ./max(m,[],1);
           imagesc(m);
        end

        
        function plotTimeCourse(tbl,conditions,neurons,fun)
            arguments
                tbl (1,1) sbx.Activity                               
                conditions 
                neurons 
                fun (1,1) = @(x)(deal(mean(x,2),std(x,0,2)./sqrt(size(x,2))));
            end
            if isa(neurons,'sbx.Roi')
                neurons = fetch(neurons);
            end

            if isa(conditions,'ns.Condition')
                conditions = fetch(conditions);
            end

            for neuron =neurons(:)'
                nexttile;
                cCntr= 0;
                thisTbl = tbl & neuron;
                T= timetable;
                for condition = conditions(:)'
                    cCntr = cCntr+1;
                    % Get the spikes Activity data for this roi
                    [y,time] = get(thisTbl & (ns.ConditionTrial & condition),'spikes');                
                    nrTime = numel(time);
                    if isempty(y);continue;end
                    y = reshape(squeeze(y),size(y,1),[]);
                
                
                    [m(:,cCntr),se(:,cCntr)] = fun(y);

                        if true
             
                
              %  plot(time,y-mean(y,1,"omitnan"),'linewidth',1,'Color',h.Color)
                        else

                    T = [T timetable(time,m,'VariableNames',{condition.name})];
                        end
                end

                m = m -mean(m,1,"omitnan");
                %m = m./se;
                grandMax = max(m,[],"all","omitnan",ComparisonMethod="abs");
                m = m./grandMax;
                se = se./grandMax;
                 m = m + repmat(1:cCntr,[nrTime 1]);
                h = ploterr(time,m,se,'linewidth',2,'ShadingAlpha',0.5);
                hold on
                hh = plot(time,repmat(1:cCntr,[nrTime 1]),'LineWidth',0.5);
                [hh.Color] =deal(h.Color);
                ylim([0 cCntr+1])
                set(gca,'yTick',1:cCntr)
                %stackedplot(T)
                 
%                 if false
% 
%                     for u=1:nrUx
%                         this = y(:,x==uX(u));
%                         this = this -mean(this,1,"omitnan");
%                         this = this./grandMax;
%                         this = this + u;
%                         plot(time,this,'LineWidth',.5,'Color',h(u).Color)
%                     end
%                     %hold on
%                     %plot(time,mean(y,2,"omitnan"),'linewidth',2)
%                     %title(sprintf('Pref = %d,
%                     %pCell=%.2f',uOri(prefOriIx),neuron.pcell))
%                 end
                %xlim(seconds([-0.25 4]))
                %
                %end
            end
        end

            function [v,time,trial,expt,roi] = get(tbl,act,pv)
                % Convenience function to extract activity (e.g., spikes) and time
                % per trial.
                % INPUT
                % tbl - sbx.Activity table
                % act - Name of a row in the sbx.ActivityParms table
                % Optional P/V pairs
                % 'time' - A datetime/duration vector with requested time
                %           points. Defaults to all in the sbx.Activity row.
                % 'trial' -  A vector of trials. Defaults to all.
                % 'method' - How to interpolate the data to the newTimes
                %                   (e.g., 'linear','mean')
                % OUTPUT
                % v - Activity per trial [nrTimePoints nrTrials]
                % t - Time in the trial (relative to the first frame of the
                %     trial) [nrTimePoints 1]
                % tr - The trials.
                % expt - The experiments
                % roi -  The rois
                arguments
                    tbl (1,1) sbx.Activity {mustHaveRows}
                    act {mustBeText}
                    pv.time (:,:) =  []
                    pv.trial (1,:) =[]
                    pv.method {mustBeText} = 'linear';
                end
                if isempty(pv.trial)
                    key = struct('act',act);
                else
                    key = struct('act',act,'trial',num2cell(pv.trial(:))');
                end
                [v,allTr,allRoi,allExpt] = fetchn(tbl & key,'activity','trial','roi','starttime');
                uTrial = unique(allTr);
                nrTrials = numel(uTrial);
                uRoi= unique(allRoi);
                nrRois= numel(uRoi);
                uExpt = unique(allExpt);
                nrExpts = numel(uExpt);
                % Rearrange by trial, expt, roi, and create ND array.
                [~,~,ix] =unique([allTr allRoi string(allExpt)],'rows','stable');
                v = reshape(v(ix),[1 nrTrials nrRois nrExpts]);
                v =cell2mat(v); % Only works if the expts had the same number of trials
                nrTimePoints =size(v,1);
                % Determine time
                parms = fetch(sbx.ActivityParms & struct('act',act),'*');
                uTime = seconds(parms.start:parms.step:parms.stop)';

                if ~isempty(pv.time)
                    allPages  = @(tf,pg)find(tf) + numel(tf)*(pg-1);
                    tt = repmat(uTime,[1 nrTrials]);
                    stay  =isbetween(tt,seconds(pv.time(:,1))',seconds( pv.time(:,2))');
                    v =reshape(v,[nrTimePoints nrTrials nrRois*nrExpts]);
                    ix = allPages(stay,1:size(v,3));
                    v = v(ix);

                    nrTimePoints =sum(stay,1);
                    trialsWithTimeInBin = nrTimePoints>0;
                    uTrial = uTrial(trialsWithTimeInBin);
                    nrTimePoints = unique(nrTimePoints(trialsWithTimeInBin));
                    assert(numel(nrTimePoints)==1,"The number of bins changes per trial...");

                    nrTrials = numel(uTrial);

                    v = reshape(v,[nrTimePoints nrTrials*nrRois*nrExpts]);
                    v = num2cell(v,1);
                    time = seconds((0:nrTimePoints-1).*parms.step)';
                    T = timetable(time,v{:});
                    T = retime(T,[min(time);max(time)],pv.method);
                    nrTimePoints = height(T);

                    v = reshape(T{:,:},[nrTimePoints nrTrials nrRois nrExpts]);
                    uTime = T.time;
                end

                time = uTime;
                trial = uTrial;
                roi = uRoi;
                expt= uExpt;

            end
        end
        methods (Access=protected)
            function makeTuples(tbl,key)
                %% Retrieve parameters for this activity
                parms = fetch(sbx.ActivityParms &key,'*');
                act = parms.act;  % Name of this parameter set
                parms = rmfield(parms,'act');
                parms = namedargs2cell(parms);
                %%  Extract the relevant timetable from sbx.ROI
                T = get(sbx.Roi & key,ns.Experiment &key,parms{:});

                %% Create tuples and insert
                [~,nrTrials] =size(T);
                tpl = struct('subject',key.subject,...
                    'session_date',key.session_date,...
                    'prep',key.prep,...
                    'roi',key.roi, ...
                    'plane',key.plane, ...
                    'starttime',key.starttime,...
                    'act',act,...
                    'trial',num2cell(1:nrTrials)',...
                    'activity',num2cell(table2array(T),1)');
                insert(tbl,tpl);
            end
        end
    end