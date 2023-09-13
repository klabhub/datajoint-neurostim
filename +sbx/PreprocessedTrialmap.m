%{
# Map triggers in a preprocessed data set to trials in Experiments and time.
->sbx.Preprocessed   # Which preprocessed data set does this apply to
->ns.Experiment
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
% numbers of frames, this is computed per Preprocessed set.
%
% This table is used by sbx.Roi to extract activity per trial, aligned to
% first frame.
%
% BK - March 2023
classdef PreprocessedTrialmap < dj.Part
    properties (SetAccess = protected)
        master = sbx.Preprocessed
    end

    methods (Access=public)
        function make(tbl,key)
            % This is called automatically from sbx.Preprocessed.makeTuples
            % (i.e. when calling populate on sbx.Preprocessed.
            expts = ns.Experiment & key;
            previousScanFrames =0;
            for expt= expts.fetch()'
                meta = fetch(ns.ExperimentMeta & expt & struct('meta_name','nrframes'),'*');
                if isempty(meta)
                    error('No nrFrames meta information found. Please sbx.addExperimentMeta.')
                end
                nrframes = meta.meta_value;
                thisFrames = previousScanFrames+ (1:nrframes);
                nrFrames= numel(thisFrames);
                previousScanFrames = thisFrames(end);

                % Read the binary file that stores TTL pulses from the
                % laser to determine laser onset time on the neurostim
                % clock.
                c = open(ns.Experiment &expt);
                [~,filename,ext] =fileparts(strrep(c.mdaq.outputFile,'\','/'));% Make fileparts os insensitive
                fldr = folder(ns.Experiment &expt);
                thisT = c.mdaq.readBin(fullfile(fldr,[filename ext]));
                laserOnIx = find(diff(thisT.laserOnDig)>0.5); % Transition from 0-1
                laserOnTime = thisT.nsTime(laserOnIx);        % Time in ns time.
                nrTTL = numel(laserOnTime);

                % Determine the time when the scanbox was instructed to
                % start/stop grabbing.
                [isGrabbing,~,~,grabbingTime]= get(c.scanbox.prms.grabbing,'withDataOnly',true);
                dtStart = seconds(grabbingTime(isGrabbing)/1000)-laserOnTime(1); %#ok<NASGU> % Time between grabbing start and first TTL
                dtStop = seconds(grabbingTime(~isGrabbing)/1000)-laserOnTime(end); % Time between grabbing stop and last TTL
                % Sanity check
                if nrFrames==nrTTL
                    % OK
                elseif nrFrames== nrTTL-1
                    % 1 trigger without a
                    % frame. Guessing it was the last.
                    laserOnTime(1)=[];
                    fprintf(2,'Removed 1 extraneous LaserOn TTL (last)\n')
                elseif nrFrames > nrTTL
                    % Some TTLs not recorded. This was a bug in the way the
                    % code first turned off the nidaq and then sent the
                    % stop command to scanbox.
                    dt = diff(laserOnTime);
                    framerate = fetch1(sbx.Preprocessed & key,'framerate');
                    SLACK = 0.05;
                    typicalDt = median(dt);
                    nrMissing = nrFrames-nrTTL;
                    if max(seconds(dt)) > (1+SLACK)/framerate
                        error('No TTL for %d frames, and frames do not appear to be successive (max dt=%s)',(nrFrames-nrTTL),max(dt));
                    elseif dtStop>0 && abs(dtStop/nrMissing-typicalDt) < SLACK*seconds(typicalDt)
                        % Nidaq  stopped too early and missed the
                        % last tiriggers. (The second clause checks that the time when grabbing was stopped
                        % on scan box matches the missing triggers within SLACK %) .
                        % Add virtual triggers at the end.
                        laserOnEstimated = laserOnTime(end) + flipud(1:nrMissing)'*typicalDt;
                        laserOnTime = [laserOnTime; laserOnEstimated];  %#ok<AGROW>
                        fprintf(2,'Apppending %d inferred laserOn TTL pulses at the end \n',nrMissing);
                    else % hasn't happened yet. Possibly prepend
                        error('dtStop = %s and the number of TTL pulses recorded by mdaq (%d) is larger than the number of frames stored by ScanBox (%d)',dtStop,nrFrames,nrTTL);
                        % laserOnEstimated = laserOnTime(1) - flipud((1:nrMissing)'*typicalDt);
                        % laserOnTime = [laserOnEstimated;laserOnTime];  %#ok<AGROW>
                        % fprintf(2,'Prepending %d\n',nrMissing);
                    end
                else
                    error('The number of TTL pulses recorded by mdaq (%d) is larger than the number of frames stored by ScanBox (%d)',nrFrames,nrTTL);
                end


                % Get the cic parameters for this experiment
                prms  = get(ns.Experiment & expt,'cic');
                if isempty(fieldnames(prms.cic))
                    % Empty cic means that first frmae was never reached.
                    % Skip this experiment
                else
                    % Split into trials
                    nrTrials = fetch1(ns.Experiment &expt,'trials');
                    % Because events are aligned to firstFrame, we do the same
                    % for the imaging data.
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
                        'tag',key.tag,...
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
end