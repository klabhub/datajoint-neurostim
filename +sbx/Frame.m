%{
# Map triggers in a session to experiments and time.
->ns.Session
->sbx.Preprocessing
-> sbx.Scan  # Which scan does this trigger belong to
trial  : int
--- 
frame       :  blob
nstime   : blob
trialtime   : blob
%}

classdef Frame < dj.Computed
    properties (Dependent)
        keySource
    end
    methods
        function v = get.keySource(~)
            % Restrict to sessions that have sbx file
            v =(ns.Session & sbx.Scan)*sbx.Preprocessing;
        end
    end

    methods (Access=protected)
        function makeTuples(tbl,key)
            scans = sbx.Scan & key;
            previousScanFrames =0;            
            for s= scans.fetch('nrframes','info')'
                thisFrames = previousScanFrames+ (1:s.nrframes);
                nrFrames= numel(thisFrames);
                previousScanFrames = thisFrames(end);

                c = open(ns.Experiment &s);
                [~,filename,ext] =fileparts(c.mdaq.outputFile);
                fldr = folder(ns.Experiment &s);
                thisT = c.mdaq.readBin(fullfile(fldr,[filename ext]));
                laserOnIx = find(diff(thisT.laserOnDig)>0.5); % Transition from 0-1
                laserOnTime = thisT.nsTime(laserOnIx);        % Time in ns time.
                
                if nrFrames== numel(laserOnTime)-1
                    % Happened in the first test set; 1 trigger without a
                    % frame. Guessing it was the last. 
                    laserOnTime(end)=[];
                end

                prms  = get(ns.Experiment & s,'cic');


                % Split into trials
                nrTrials = fetch1(ns.Experiment &s,'trials');
                % Becuase events are aligned to firstFrame, we do the same
                % for the spiking data.
                trialStartTime = seconds(prms.cic.firstFrameNsTime/1000);
                stay = cell(1,nrTrials);
                PRE =seconds(0.5);
                for tr=1:nrTrials
                    stay{tr} = laserOnTime >=(trialStartTime(tr)-PRE);
                    if tr<nrTrials
                        stay{tr} = stay{tr} & laserOnTime <(trialStartTime(tr+1)-PRE);
                    end
                end
                frames = cellfun(@(x)(thisFrames(x)),stay,'uni',false)';
                clockTimes = cellfun(@(x)(seconds(laserOnTime(x))*1000),stay,'uni',false)';                
                trialTimes= cellfun(@(x,y) (seconds(laserOnTime(x)-y)*1000),stay,num2cell(trialStartTime)','uni',false)';



                tpl = struct('subject',key.subject,...
                    'session_date',key.session_date,...
                    'name',key.name,...
                    'starttime',s.starttime,...
                    'trial',num2cell(1:nrTrials)',...
                    'frame',frames, ...
                    'nstime',clockTimes,...
                    'trialtime',trialTimes);
                insert(tbl,tpl)
            end
        end
    end
end