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
        function [v,time,tr,expt,roi] = get(tbl,act,pv)
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
                tbl (1,1) sbx.Activity
                act {mustBeText}
                pv.time (1,:) =  []
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
             [uRows,~,ix] =unique([allTr allRoi string(allExpt)],'rows','stable');
             v = reshape(v(ix),[1 nrTrials nrRois nrExpts]);

             v =cell2mat(v); % Only works if the expts had the same number of trials
             % Determine time 
             parms = fetch(sbx.ActivityParms & struct('act',act),'*');
             time = seconds(parms.start:parms.step:parms.stop)';
            nrTimePoints= numel(time);
             tr = uTrial;
             if ~isempty(pv.time)
                 v =reshape(v,nrTimePoints,[]);
                 v =num2cell(v,1);
                 T = timetable(time,v{:});
                 T = retime(T,seconds(pv.time(:)),pv.method);
                 time = T.time;
                 nrTimePoints =numel(time);
                 v = reshape(T{:,:},[nrTimePoints nrTrials nrRois nrExpts]);
             end
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