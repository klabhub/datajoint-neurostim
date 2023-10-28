%{
# Map frames in a movie to trials and time in an Experiment
-> ns.Movie     # The Movie this applies to
trial  :int     # The trial number
--- 
frame    = NULL    : blob   # Frames from the movie that correspond to this trial
nstime   = NULL    : blob   # Time in seconds on the Neurostim experiment clock
trialtime = NULL   : blob   # Time in seconds relatve to the first frame in the trial.
%}
%
% A frame is assigned to a trial if it occurs at or after
% the start of the first frame in the trial and before the first frame of the next trial.
% (In other words the ITI is included at the end of each trial).
%
% BK - Oct 2023
classdef MovieTrialmap < dj.Part
    properties (SetAccess = protected)
        master = ns.Movie
    end

    methods (Access=public)
        function make(tbl,key)
            % This is called automatically from ns.Movie.makeTuples
            % (i.e. when calling populate on ns.Movie).

            % To link with trials and times in an experiment we need to determine the time of
            % the first frame in nsTime. If NS recorded the movie this can be extracted from the plugin parameters,
            % but I haven;t implemented that yet as we  do not have such
            % movies at this time. (Only those recorded by SBX)
            
            mvTpl = fetch(ns.Movie&key,'*');
            if contains(key.filename,{'_eye' ,'_ball'})                
                % An sbx eye tracker or ball tracker movie. Movie frames are
                % supposed to match the frames of the TPI.
                tpi = sbx.Preprocessed & (ns.Session & key);
                if exists(tpi)
                    tpi =fetch(tpi,'framerate','LIMIT 1');% Assuming that every preprocessing has the same trialmap.
                    tpiTpl = fetch(sbx.PreprocessedTrialmap &tpi &key,'*'); % Get all trial maps
                    if isempty(tpiTpl)
                        fprintf(2,"No PrepocessedTrialmap for this experiment (%s on %s@%s). Too few trials? (%d)\n",mvTpl.subject,mvTpl.session_date,mvTpl.starttime,fetch1(ns.Experiment & mvTpl,'trials'));
                        return;
                    end
                    nrFramesPerTrial2PI = cellfun(@numel,{tpiTpl.frame});
                    SLACK = 0.01;% Allow this much slack (as a fraction of the frame count to allow for a few extra frames in the movie)
                    if abs(sum(nrFramesPerTrial2PI)-mvTpl.nrframes)/mvTpl.nrframes<SLACK
                        % Match; reuse the trialmap for the tpi. But the
                        % TPI trialmap is across the session, while movies
                        % are per experiment. So we subtract the (session-based) frame
                        % number of the first trial in this experiment
                        offset = tpiTpl(1).frame(1)-1;
                        if offset>0
                            for i=1:numel(tpiTpl)
                                tpiTpl(i).frame = tpiTpl(i).frame-offset;
                            end
                        end
                        mvTrialMapTpl = rmfield(tpiTpl,{'tag'});
                        [mvTrialMapTpl.filename]  =deal(mvTpl.filename);
                        insert(tbl,mvTrialMapTpl);
                        % Update the mvTpl to match the framerate of the TPI
                        update(ns.Movie &mvTpl,'framerate',tpi.framerate);
                        fprintf('Updating %s with framerate matching TPI acquisition.\n',mvTpl.filename)
                        return;  % Done
                    else
                        % Different number of frames. I think this happens when
                        % framesPerTrigger is not set to 1 in SBX. Assume the
                        % first frame was triggered the same as the TPI, and
                        % that all other frames followed at the mv framerate.
                        % (The alternative is that frames are triggered by
                        % the TPI but that two are collected each time (beacuse
                        % the camera was set to 30 Hz. Because the TPI is at
                        % 15.5Hz and not exactly 15 this could result in
                        % slightly different timing).
                        firstTrial = tpiTpl([tpiTpl.trial]==1);
                        firstFrameNsTime = firstTrial.nstime(1);
                    end
                else
                    % No entry yet in the sbx.Preprocessed table. Generate
                    % an error so that a future populate call on ns.Movie
                    % will retry this.
                    error('For SBX _ball and _eye movies the Preprocessed table needs to be populated first. Rerun later.\n')                    
                end
            elseif exist(ns.Plugin & key & 'plugin_name=''camera''')
                % A movie captured by the neurostim camera plugin
                   prms  = get(ns.Experiment & key,'camera','prm','firstVideoFrame','atTrialTime',inf,'what');
            if isempty(fieldnames(prms.cic))
                % Empty cic means that first frame was never reached.
                % Skip this experiment
            else
                firstFrameNsTime = prms.cic.firstFrameNsTime(1);
            end

            end

            % Assume regular sampling at the framerate.
            frameTime = seconds((0:mvTpl.nrframes-1)/mvTpl.framerate+firstFrameNsTime);

            % Split into trials
            nrTrials = fetch1(ns.Experiment &key,'trials');
            % Because events are aligned to firstFrame, we do the same
            % for the movie frames.
            trialStartTime = seconds(prms.cic.firstFrameNsTime/1000);
            stay = cell(1,nrTrials);
            for tr=1:nrTrials
                stay{tr} = frameTime >=trialStartTime(tr);
                if tr<nrTrials
                    stay{tr} = stay{tr} & frameTime < trialStartTime(tr+1);
                end
            end
            thisFrames= 1:mvTpl.nrframes;
            frames = cellfun(@(x)(thisFrames(x)),stay,'uni',false)';
            nsTimes = cellfun(@(x)(seconds(frameTime(x))),stay,'uni',false)';
            trialTimes= cellfun(@(x,y) (seconds(frameTime(x)-y)),stay,num2cell(trialStartTime)','uni',false)';
            nrPerTrial = cellfun(@numel,frames);
            if any(nrPerTrial==0)
                fprintf(2,'No movie frames for %d trials, %.0f on average in %s.\n',sum(nrPerTrial==0),mean(nrPerTrial),mvTpl.filename);
            end
            % Create tuples and insert.
            mvKey = namedargs2cell(ns.stripToPrimary(tbl,key));
            partTpl = struct(mvKey{:},'trial',num2cell(1:nrTrials)',...
                'frame',frames, ...
                'nstime',nsTimes,...
                'trialtime',trialTimes);
            insert(tbl,partTpl)
        end
    end
end
