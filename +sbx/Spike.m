%{
# The (deconvolved) spiking activity of an ROI in one trial of an experiment
-> sbx.Scan
trial :smallint
-> sbx.Roi
---
activity : longblob   # Deconvolved spiking activity trace
%}

classdef Spike < dj.Computed

     properties (Dependent)
        keySource
    end
    methods
        function v = get.keySource(~)
            % Restrict to experiments that have an sbx file
            v = sbx.Scan;
        end
    end
    methods (Access=protected)
        function makeTuples(tbl,key)
            %% Find scanner trigger just before trial start event

            
            nrTrials = fetch1(ns.Experiment &key,'trials');
            prms  = get(ns.Experiment & key,'cic');
            % Approximate solution for first m60 trial
            trialStartTime = prms.cic.trialClockTime(prms.cic.trial>0);
            frameInfo = fetch(sbx.Frame &key,'*');
            stay = cell(1,nrTrials);
            for tr=1:nrTrials
                stay{tr} = frameInfo.clocktime >=trialStartTime(tr);
                if tr<nrTrials
                    stay{tr} = stay{tr} & frameInfo.clocktime <trialStartTime(tr+1);
                end
            end
            

            % Fetch all spikes for this roi in this session
            roiCntr= 0;
            nrRois = count(sbx.Roi & key);
            spikes= cell(nrRois,nrTrials);
            %times = cell()
            rois=  fetch(sbx.Roi & key,'spikes');
            for roi = rois'            
                roiCntr = roiCntr+1;
                [spikes(roiCntr,:)] = cellfun(@(x)(roi.spikes(frameInfo.frame(x))),stay,'uni',false);
               % time{tr} = frameInfo.clocktime(stay)'; % ACtual time of triggers.
            end

            
            for tr=1:nrTrials
            
            tpl = struct('subject',key.subject,...
                            'session_date',key.session_date,...
                            'name',{rois.name}',...
                            'roi',{rois.roi}', ...
                            'plane',{rois.plane}', ...
                            'starttime',key.starttime,...
                            'trial',tr,...
                            'activity',spikes(:,tr));          
            tr
                    insert(tbl,tpl);
            end
         

        end
    end
end