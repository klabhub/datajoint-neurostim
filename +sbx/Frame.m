%{
# Map triggers in a session to trials in Experiments and time.
->sbx.Preprocessed   # Which preprocessed data set does this apply to
->sbx.ExperimentSbx  # Which scan does this trigger belong to
trial  : int
--- 
frame       :  blob   # Frames from this Preprocessed set that correspond to this trial
nstime   : blob       # Time in seconds on the Neurostim experiment clock
trialtime   : blob    # Time in seconds relatve to the first frame in the trial.
%}
% This determines how each frame in a session maps to a trial and a time in
% a Neurostim experiment. A frame is assigned to a trial if it occurs less
% than 500 ms before the first frame in the trial and more than 500 ms
% before the first frame of the next trial. (This 500 ms is hardcoded as PRE).
% 
% Because different sbx.Preprocessed data sets for the same ns.Session could have different
% numbers of frames, this is computed per Preprocessed set. (Although no
% functionality currently exists to preprocess only a subset of
% experiments, except by deleting the files).
% 
% This table is used by sbx.Roi to extract activity per trial, aligned to
% first frame. 
% BK - March 2023
classdef Frame < dj.Imported

    methods (Access=protected)
        function makeTuples(tbl,key)
            PRE =seconds(0.5);  % Keeping 500 ms before first frame with each trial
            expts = sbx.ExperimentSbx & key;
            previousScanFrames =0;            
            for expt= expts.fetch('nrframes','info')'
                thisFrames = previousScanFrames+ (1:expt.nrframes);
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
                else
                    error('THe number of TTL pulses recorded by mdaq (%d) does not match the frames stored by ScanBox (%d)',nrFrames,nrTTL);
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
                    stay{tr} = laserOnTime >=(trialStartTime(tr)-PRE);
                    if tr<nrTrials
                        stay{tr} = stay{tr} & laserOnTime <(trialStartTime(tr+1)-PRE);
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