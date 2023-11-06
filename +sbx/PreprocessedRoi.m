%{
# Properties of a ROI in a session, based on a Preprocessed set.
-> sbx.Preprocessed
roi : smallint    #  id within this segmentation and plane
plane : smallint #  plane
---
pcell = 0 : float # Probability that the ROI is a neuron
x        : Decimal(4,0) # x-Pixel location 
y        : Decimal(4,0) # y-Pixel location 
radius   : Decimal(4,1) # Radius of cell in micron
aspect   : Decimal(4,1)  # Aspect raio of major/minor axis of a 2D Gaussian fit.
compact  : Decimal(4,2)  # How compact the ROI is ( 1 is a disk, >1 means less compact)
%}
classdef PreprocessedRoi < dj.Part
     properties (SetAccess = protected)
        master = sbx.Preprocessed
    end


    methods (Access=?sbx.Preprocessed)
        function make(tbl,key,fldr,micPerPix)
            % Read the npy results from the suit2p folder and store them in
            % the table.

            planes = dir(fullfile(fldr,'plane*'));
           
            %% Get calibration info
         
            for pl = 1:numel(planes)
                %% Read npy           
               thisFile = fullfile(fldr,planes(pl).name,'iscell.npy');
               if ~exist(thisFile,"file")
                        error('File %s does not exist',thisFile);
               end               
               iscell = single(py.numpy.load(thisFile,allow_pickle=true));                
                
               % We saved the stat.npy as stat.mat in sbx.Preprocessed
                thisFile = fullfile(fldr,planes(pl).name,'stat.mat');
               if ~exist(thisFile,"file")
                        error('File %s does not exist',thisFile);
                end
                load(thisFile,'stat');
                stat= [stat{:}];
                med= cat(1,stat.med); %[y x] center pixels per ROI.
                compact = cat(1,stat.compact);
                aspect = cat(1,stat.aspect_ratio);

                radius = cat(1,stat.radius); % Pixels
                radius = radius.*micPerPix;
                
                %% Make tuples and insert=
                nrROIs= numel(aspect); 
                
                tpl = mergestruct(key, ...
                    struct('roi',num2cell(1:nrROIs)', ...
                    'plane',pl-1, ...
                    'pcell',num2cell(iscell(:,2)), ...
                    'x',num2cell(med(:,2)), ...    % Dim 2 is the horizontal axis of the image
                    'y',num2cell(med(:,1)), ...    % Dim 1 is the vertical axis of the image
                    'radius',num2cell(radius',1)', ...
                    'compact',num2cell(compact',1)', ...
                    'aspect',num2cell(aspect',1)'));                
                insert(tbl,tpl);                
            end
        end
    end
end