%{
# Properties of a ROI in a session, based on a Preprocessed set.
-> ns.CChannel
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
classdef Roi < dj.Computed
    properties (Dependent)
        keySource
    end
    methods
        function v = get.keySource(~)
            % The key source is the ns.C table becuase we read all channels
            % for a ns.C entry from a single npy file 
            v = ns.C & (ns.CParm & 'extension=''.sbx''');
        end
    end
     
    methods (Access=protected)
        function makeTuples(tbl,cKey)            
            % Extract info from Channel key
            c = fetch((ns.C & cKey)*(ns.CParm & 'extension=''.sbx'''),'*');
            if isempty(c)
                error('This ns.C (%s) does not contain SBX data. Cannot create sbx.ROI tuples.',cKey.ctag)
            end
            prep = sbx.Preprocessed & struct('subject',c.subject,'session_date',c.session_date,'prep',c.parms.prep);
            info = sbx.readInfoFile(fetch(ns.Experiment&cKey,'LIMIT 1'));
            micPerPix = sqrt(sum([info.xscale info.yscale].^2));
            if count(prep)~=1
                error('Need exactly 1 sbx.Preprocessed set to extract ROIs, not %d',count(prep))
            end
            fldr = fullfile(folder(ns.Experiment & cKey),fetch1(prep,'folder'));
            if ~exist(fldr,"dir")
                error('Preprocessed data folder %s not found',fldr)
            end
            planes = dir(fullfile(fldr,'plane*'));
            channelsSoFar = 0;       
            
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
                key = repmat(ns.stripToPrimary(ns.C,cKey),[nrROIs 1]);
                channelsPerRoi = num2cell(channelsSoFar+(1:nrROIs));
                [key.channel] = deal(channelsPerRoi{:});
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
                channelsSoFar = channelsSoFar + nrROIs; % Next plane starts at nrROIs+1
            end
        end
    end
end