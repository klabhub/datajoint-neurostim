%{
# Movie - A table listing movies associated with an experiment
-> ns.File
---
nrframes=null               : int unsigned                  # Number of frames in the video
framerate=null              : float                         # Framerate (fps)
width=null                  : int unsigned                  # Width in pixelsststr
height=null                 : int unsigned                  # Height in pixels
bytes=null                  : double                        # Number of bytes in the file
intensity=null              : longblob                # Mean intensity per frame
%}

classdef Movie < dj.Computed
    properties (Constant)
        extensions = {".mj2" ,".avi" ,".mp4" ,".mpeg4" ,".mov"};
    end
    properties (Dependent)
        keySource
    end

    methods
        function v = get.keySource(~)
            v =ns.File & struct('extension',ns.Movie.extensions) ;
        end
    end 
    methods (Access=public)

        function movie=open(tbl,sizeOption)
            % Open movies as videoreader objects
            arguments 
                tbl (1,1) ns.Movie
                sizeOption  (1,1) string {mustBeMember(sizeOption,["smallest","largest","all"])} ="all"
            end           
            mvFile =  file(tbl,sizeOption);
            if ischar(mvFile) || isstring(mvFile)
                movie = VideoReader(mvFile); %single movie
            else
                for i=1:numel(mvFile)
                    movie{i} = VideoReader(mvFile{i}); %#ok<TNMLP,AGROW>
                end
            end
        end


        function [v,bytes]= file(tbl,sizeOption)
            % Return the full path to the movie. It is possible that there is more than one file; pick the smallest one
            % (sizeOption =-1), all files (0), or the largest one (+1).(Thi
            % sis used, for instance, when a compressed and non-compressed
            % version of a file exist). 
            arguments 
                tbl (1,1) ns.Movie
                sizeOption  (1,1) string  {mustBeMember(sizeOption,["smallest","largest","all"])} ="all"
            end
            % Get the list
            tpl = fetch(tbl,'bytes','ORDER BY bytes');    
            nrF = numel(tpl);
            % Pick the requested ones
            switch sizeOption
                case "smallest"
                    ix =1;
                case "largest"
                    ix = nrF;
                case "all"
                    ix = 1:nrF;
            end
            v = strings(numel(ix),1);
            iCntr= 0;
            for i=ix
                iCntr =iCntr+1;
                v(iCntr) = fullfile(folder(ns.Experiment &tpl(i)),tpl(i).filename);
            end
            if nargout >1
                bytes= [tpl(ix).bytes]';
            end
        end
    end

    methods (Access=protected)
        function makeTuples(tbl,key)
            fldr = folder(ns.Experiment& key);
            ff =fullfile(fldr,key.filename);
            if exist(ff,'file')
                info =dir(ff);                
                try
                    fprintf('Opening %s \n',ff)
                    mv= VideoReader(ff); 
                catch me
                    fprintf('Could not open %s.\n Error: %s\n', ff,me.message)
                    return;
                end
            else
                error('File not found %s',ff);
            end
            fprintf('Reading %d frames to determine mean intensity...',mv.NumFrames)
            frames= read(mv);            
            fprintf('Done.\n')
            intensity = squeeze(mean(frames,1:(ndims(frames)-1),"omitnan"));
            tpl = mergestruct(key,struct('nrframes',mv.NumFrames,'width',mv.Width,'height',mv.Height,'framerate',mv.FrameRate,'bytes',info.bytes,'intensity',intensity));        
            insert(tbl,tpl)

            % Now fill the MovieTrialmap
            make(ns.MovieTrialmap,key);
        end
    end

end
