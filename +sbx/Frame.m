%{
# Map triggers in a session to trials in Experiments and time.
->sbx.Preprocessed   # Which preprocessed data set does this apply to
->ns.Session         # Which session does this trigger belong to
->ns.Experiment      # Which experiment does this trigger belong to
trial  : int
--- 
frame       :  blob   # Frames from this Preprocessed set that correspond to this trial
nstime   : blob       # Time in seconds on the Neurostim experiment clock
trialtime   : blob    # Time in seconds relatve to the first frame in the trial.
%}
% This determines how each frame in a session maps to a trial and a time in
% a Neurostim experiment. A frame is assigned to a trial if it occurs after
% the first frame in the trial and before the first frame of the next trial. 
% (In other words the ITI is included at the end of each trial).
% 
% Because different sbx.Preprocessed data sets for the same ns.Session could have different
% numbers of frames, this is computed per Preprocessed set. (Although no
% functionality currently exists to preprocess only a subset of
% experiments, except by deleting the experiment files).
% 
% This table is used by sbx.Roi to extract activity per trial, aligned to
% first frame. 
% BK - March 2023
classdef Frame < dj.Imported
     properties (Dependent)
         keySource
     end
 
     methods
         function v= get.keySource(tbl)
             v =   ns.Session * sbx.Preprocessed; 
         end
     end

    methods (Access=protected)
        function makeTuples(tbl,key)            
            expts = ns.Experiment & key;
            previousScanFrames =0;            
            for expt= expts.fetch()'
                meta = fetch(ns.ExperimentMeta & expt & struct('meta_name','nrframes'),'*');
                nrframes = meta.meta_value;                
                thisFrames = previousScanFrames+ (1:nrframes);
                nrFrames= numel(thisFrames);
                previousScanFrames = thisFrames(end);
                % Read the binary file that stores TTL pulses from the
                % laser to determine laser onset time on the neurostim
                % clock.
                c = open(ns.Experiment &expt);
                [~,filename,ext] =fileparts(c.mdaq.outputFile);
                fldr = folder(ns.Experiment &expt);
                thisT = c.mdaq.readBin(fullfile(fldr,[filename ext]));
                laserOnIx = find(diff(thisT.laserOnDig)>0.5); % Transition from 0-1
                laserOnTime = thisT.nsTime(laserOnIx);        % Time in ns time.
                nrTTL = numel(laserOnTime);
                % Sanity check
                if nrFrames==nrTTL
                    % OK
                elseif nrFrames== nrTTL-1
                    % Happened in the first test set; 1 trigger without a
                    % frame. Guessing it was the last. 
                    laserOnTime(end)=[];
                    fprintf(2,'Removed 1 extraneous LaserOn TTL (last)\n')
                elseif nrFrames > nrTTL
                    dt = diff(laserOnTime);
                    sd  = std(dt);
                    typicalDt = median(dt);
                    nrMissing = nrFrames-nrTTL;
                    framerate = fetch1(sbx.Preprocessed & key,'framerate');
                    SLACK = 0.05;
                    if max(seconds(dt)) > (1+SLACK)/framerate
                        error('No TTL for %d frames, and frames do not appear to be successive (max dt=%s)',(nrFrames-nrTTL),max(dt));
                    else
                        fprintf(2,'Trying to fix %d missing TTL by reading eye data\n',nrFrames-nrTTL)
                        % The frames for which we have recorded TTLs appear to be consecutive
                        % (within SLACK *framerate), so the missing TTL are likely either at the start 
                        % or the end. 
                      
                        ttlStartTime = thisT.clockTime(1);
                        ttlStopTime = thisT.clockTime(end);
                        % The clock on the scanbox computer can be read
                        % from the _eye file (or _ball, but _eye is
                        % smaller). I tried speeding this up with partial
                        % loading using matfile, but that made things even
                        % slower (probably the .mat is an older style)
                        mf = load(fullfile(fldr,filename,[filename '__001_eye.mat']));
                        scanboxStartTime = datetime(mf.abstime(1,1).AbsTime);
                        scanboxStopTime = datetime(mf.abstime(end,1).AbsTime); % First and last
                        % Compare the clocks (This assumes that the
                        % Neurostim and the Scanbox computer have a clock
                        % that is at least approximately synchronized (e.g., in the OS).
                        offsetAtStart = abs(scanboxStartTime-ttlStartTime) ;
                        offsetAtStop  = abs(scanboxStopTime-ttlStopTime);                         
                        if offsetAtStop < offsetAtStart
                            % The stop times match best. Prepend.
                            laserOnEstimated = laserOnTime(1) - flipud(cumsum(1:nrMissing)'*typicalDt);
                            laserOnTime = [laserOnEstimated;laserOnTime];  %#ok<AGROW> 
                            fprintf(2,'TTL acquisition probably started too late (%s<%s). Prepending %d estimated laser onset times (Laser SD = %.2f ms)\n',offsetAtStop,offsetAtStart,nrMissing,seconds(sd)*1000)
                        else 
                            % The start times match best. Append.
                            laserOnEstimated = laserOnTime(end) + cumsum(1:nrMissing)'*typicalDt;
                            laserOnTime = [laserOnTime; laserOnEstimated];  %#ok<AGROW> 
                            fprintf(2,'TTL acquisition probably  stopped too early (%s<%s). Appending %d estimated laser onset times (Laser SD = %.2f ms)\n',offsetAtStart,offsetAtStop,nrMissing,seconds(sd)*1000)
                        end   
                        % Diagnostics
%                         [grabbing,~,~,sbxGrabTime] = get(c.scanbox.prms.grabbing,'withDataOnly',true);
%                         grabStartTime= seconds(sbxGrabTime(grabbing)/1000);
%                         grabStopTime = seconds(sbxGrabTime(~grabbing)/1000);
%                         nrExpected = round(framerate*seconds(grabStopTime-grabStartTime));

                    end
                else
                    error('The number of TTL pulses recorded by mdaq (%d) is larger than the number of frames stored by ScanBox (%d)',nrFrames,nrTTL);
                end
                % Get the cic parameters for this experiment
                prms  = get(ns.Experiment & expt,'cic');
                % Split into trials
                nrTrials = fetch1(ns.Experiment &expt,'trials');
                % Because events are aligned to firstFrame, we do the same
                % for the spiking data.
                trialStartTime = seconds(prms.cic.firstFrameNsTime/1000);
                stay = cell(1,nrTrials);               
                for tr=1:nrTrials
                    stay{tr} = laserOnTime >=trialStartTime(tr);
                    if tr<nrTrials
                        stay{tr} = stay{tr} & laserOnTime < trialStartTime(tr+1);
                    end
                end
                frames = cellfun(@(x)(thisFrames(x)),stay,'uni',false)';
                nsTimes = cellfun(@(x)(seconds(laserOnTime(x))),stay,'uni',false)';                
                trialTimes= cellfun(@(x,y) (seconds(laserOnTime(x)-y)),stay,num2cell(trialStartTime)','uni',false)';
                % Create tuples and insert.
                tpl = struct('subject',key.subject,...
                    'session_date',key.session_date,...
                    'prep',key.prep,...
                    'starttime',expt.starttime,...
                    'trial',num2cell(1:nrTrials)',...
                    'frame',frames, ...
                    'nstime',nsTimes,...
                    'trialtime',trialTimes);
                insert(tbl,tpl)
            end
        end
    end
end