%{
# ROI in a session with is Fluorescence and spiking activity.
-> sbx.Segmentation
roi : smallint    #  id within this segmentation and plane
plane : smallint #  plane
---
pcell =0 : float # Probability that the ROI is a neuron
fluorescence   : longblob     # Fluorescence trace
neuropil: longblob   # Neuropil Fluorescence trace (scatter)            
spikes: longblob   # Deconvolved spiking activity
%}
classdef Roi < dj.Computed

    methods (Access = public)

        function [varargout] = get(tbl,expt,pv)
            arguments
                tbl (1,1) sbx.Roi
                expt (1,1) sbx.Scan
                pv.modality = 'spikes'
                pv.trial = []
                pv.start =-500
                pv.stop  =2500
                pv.interpolation = 'nearest';
                pv.step = 100;
            end

            frame = sbx.Frame & tbl & expt;
            if ~isempty(pv.trial)
                frame = frame & struct('trial',num2cell(pv.trial)');
            end
            sessionV= fetch(tbl,pv.modality); % Values (e.g., spikes) across session
            V = [sessionV.(pv.modality)];
            V = mean(V,2,"omitnan");
            frames= fetch(frame,'*');
            trialCntr =0;
            nrTrials = count(frame);
            newTimes = seconds((pv.start:pv.step:pv.stop)/1000);
            T = timetable;
            for f = frames'
                trialCntr = trialCntr +1;    
                
                thisT = timetable(seconds(f.trialtime/1000),V(f.frame,:),'VariableNames',"Trial" + string(trialCntr));                
                inWindow  =thisT.Time>=seconds(pv.start/1000) & thisT.Time<seconds(pv.stop/1000);
                if trialCntr ==1
                    T= thisT(inWindow,:);
                else
                    T = synchronize(T,thisT(inWindow,:) , ...
                                    newTimes, ...
                                    pv.interpolation, ...                                    
                                    'EndValues',NaN);
                end
            end
           
                      
        if nargout ==2
            varargout{1} = seconds(T.Time)*1000;
            varargout{2} = double(T{:,1:end});
        else
            varargout{1} =T;
        end


        end
    end

    methods (Access=protected)
        function makeTuples(tbl,key)
            fldr= getFolder(sbx.Segmentation & key);            
            planes = dir(fullfile(fldr,'plane*'));
            for pl = 1:numel(planes)
                files= {'iscell','F','Fneu','spks'};
                vals=  cell(1,numel(files));
                fprintf('Reading numpy files...\n')
                for f=1:numel(files)
                    thisFile = fullfile(fldr,planes(pl).name,[files{f} '.npy']);
                    if ~exist(thisFile,"file")
                        error('File %s does not exist',thisFile);
                    end
                    vals{f}= single(py.numpy.load(thisFile,allow_pickle=true));
                end
                fprintf('Done.\n')
                
                [iscell,f,fneu,spks] = deal(vals{:});

                nrROIs= size(iscell,1);                
                tpl = struct('subject',key.subject,...
                            'session_date',key.session_date,...
                            'name',key.name,...
                            'roi',num2cell(1:nrROIs)', ...
                            'plane',pl-1, ...
                            'pcell',num2cell(iscell(:,2)), ...
                            'fluorescence',num2cell(f',1)',...
                            'neuropil',num2cell(fneu',1)',...
                            'spikes',num2cell(spks',1)');
                CHUNK =350;
                for i=1:CHUNK:numel(tpl)
                    i
                    insert(tbl,tpl(i:min(i+CHUNK-1,numel(tpl))));
                end
            end
        end
    end
end