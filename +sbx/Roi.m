%{
# ROI in a session with is Fluorescence and spiking activity.
-> sbx.Preprocessed
roi : smallint    #  id within this segmentation and plane
plane : smallint #  plane
---
pcell =0 : float # Probability that the ROI is a neuron
fluorescence   : longblob     # Fluorescence trace
neuropil: longblob   # Neuropil Fluorescence trace (scatter)            
spikes: longblob   # Deconvolved spiking activity
meanrate : float  # Mean number of spikes/second across the session
stdrate  :  float  # Stdev of the spikes/second across the session
x        : Decimal(4,0) # x-Pixel location 
y        : Decimal(4,0) # y-Pixel location 
%}
classdef Roi < dj.Imported

    methods (Access = public)

        function [varargout] = get(roi,expt,pv)
            % Function to retrieve trial-start aligned activity data per
            % ROI and Experiment. 
            % roi  - Table o sbx.Roi to use
            % expt - Table of sbx.ExperimentSbx to use
            % Optional Parameter/Value pairs
            % modality - 'spikes','fluorescence','neuropil'
            % trial   - Which trials to extract
            % start  - First time point to extract (relatve to first frame
            %           of each trial, in seconds)
            % stop   - Last time point to extract 
            % step   - Step size in seconds.
            % interpolation -  enum('nearest','linear','spline','pchip','makima')  
            %               Interpolation method; see timetable/synchronize. ['linear']               
            % OUTPUT
            %  [t,v]  = t: time in seconds since first frame event, 
            %           v: Matrix with time along the rows, and trials
            %           along columns.
            % Alternatively, when only a single output is requested:
            % T     = timetable with each column a trial. Time is in seconds 
            %           relative to the first frame of the trial.
            arguments
                roi (1,1) sbx.Roi
                expt (1,1) sbx.ExperimentSbx
                pv.modality = 'spikes'
                pv.trial = []
                pv.start =-0.5
                pv.stop  =2.5
                pv.interpolation {mustBeText}%, mustBeMember(pv.interpolation,{'linear','nearest','spline','pchip','makima'})}= 'linear';
                pv.step = 0.1;
            end
            
            %% Get the mapping from Frames to trials.
            % Specific or this ROI (i.e. this Preprocessed set) in this
            % Expt
            frame = sbx.Frame & roi & expt;
            if ~isempty(pv.trial)
                frame = frame & struct('trial',num2cell(pv.trial)');
            end
            frames= fetch(frame,'*');

            % Retrieve the activity in the entire session 
            sessionActivity= fetch(roi,pv.modality); % Values (e.g., spikes) across session
            V = [sessionActivity.(pv.modality)]; %[nrFramesPerSession nrROIs]
            
            % Define the new time axis (time relative to firstFrame event).
            newTimes = seconds(pv.start:pv.step:pv.stop);            
            trialCntr =0;
            for f = frames'
                trialCntr = trialCntr +1;                    
                thisT = timetable(seconds(f.trialtime),V(f.frame,:),'VariableNames',"Trial" + string(trialCntr));                
                 if trialCntr ==1
                    T= retime(thisT,newTimes,pv.interpolation,'EndValues',NaN);
                else
                    T = synchronize(T,thisT, newTimes, pv.interpolation, 'EndValues',NaN);
                end
            end
           
                      
        if nargout ==2
            varargout{1} = seconds(T.Time);
            varargout{2} = double(T{:,1:end});
        else
            varargout{1} =T;
        end


        end
    end

    methods (Access=protected)
        function makeTuples(tbl,key)
            % Read the npy results from the suit2p folder and store them in
            % the table.            
            CHUNK =350;  % This many ROIs are sent to the server at the same time.                          

            fldr= getFolder(sbx.Preprocessed & key);            
            planes = dir(fullfile(fldr,'plane*'));
            files= {'iscell','F','Fneu','spks'};                
            for pl = 1:numel(planes)
                %% Read npy
                vals=  cell(1,numel(files));
                fprintf('Reading numpy files...\n')
                for f=1:numel(files)
                    thisFile = fullfile(fldr,planes(pl).name,[files{f} '.npy']);
                    if ~exist(thisFile,"file")
                        error('File %s does not exist',thisFile);
                    end
                    vals{f}= single(py.numpy.load(thisFile,allow_pickle=true));
                end
                thisFile = fullfile(fldr,planes(pl).name,'ops.npy');
                ops = py.numpy.load(thisFile,allow_pickle=true);
                % We saved the stat.npy as stat.mat in sbx.Preprocessed
                thisFile = fullfile(fldr,planes(pl).name,'stat.mat');
                load(thisFile,'stat');
                stat= [stat{:}];
                med= cat(1,stat.med); %[x y] pixels per ROI.               
                fprintf('Done.\n')
                %% Make tuples and insert
                [iscell,f,fneu,spks] = deal(vals{:});
                
                [nrROIs,nrFrames] = size(f); %#ok<ASGLU> 
                frameDuration = 1./ops.item{'fs'};
                meanrate = mean(spks,2,"omitnan")/frameDuration;
                stdrate = std(spks,0,2,"omitnan")/frameDuration;
                tpl = struct('subject',key.subject,...
                            'session_date',key.session_date,...
                            'prep',key.prep,...
                            'roi',num2cell(1:nrROIs)', ...
                            'plane',pl-1, ...
                            'pcell',num2cell(iscell(:,2)), ...
                            'fluorescence',num2cell(f',1)',...
                            'neuropil',num2cell(fneu',1)',...
                            'spikes',num2cell(spks',1)', ...
                            'meanrate',num2cell(meanrate',1)',...
                             'stdrate', num2cell(stdrate',1)', ...
                             'x',num2cell(med(:,1)), ...
                             'y',num2cell(med(:,2)));                               
                fprintf('Adding %d ROIs to the database in chunks of %d\n',nrROIs,CHUNK);                
                tic;
                for i=1:CHUNK:numel(tpl)                  
                    insert(tbl,tpl(i:min(i+CHUNK-1,numel(tpl))));
                end
                fprintf('Done in %s \n',seconds(toc))
            end
        end
    end
end