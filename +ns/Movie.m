%{
# Movie - A table listing movies associated with an experiment
-> ns.File 
---
nrframes = NULL : int unsigned      # Number of frames in the video
framerate = NULL : float            # Framerate (fps)
width = NULL : int unsigned         # Width in pixelsststr
height = NULL : int unsigned        # Height in pixels
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
                sizeOption  (1,1) double {mustBeMember(sizeOption,[-1 0 1])} =0
            end           
            mvFile =  file(tbl,sizeOption);
            if ischar(mvFile);mvFile={mvFile};end
            for i=1:numel(mvFile)
                movie(i) = VideoReader(mvFile{i}); %#ok<TNMLP,AGROW>
            end
        end


        function v= file(tbl,sizeOption)
            % Return the full path to the movie. It is possible that there is more than one file; pick the smallest one
            % (sizeOption =-1), all files (0), or the largest one (+1).(Thi
            % sis used, for instance, when a compressed and non-compressed
            % version of a file exist). 
            arguments 
                tbl (1,1) ns.Movie
                sizeOption  (1,1) double {mustBeMember(sizeOption,[-1 0 1])} =0
            end
            % Determine existence and size
            files = cell(count(tbl));
            bytes = nan(count(tbl));
            fCntr=0;
            for f=fetch(tbl,'filename')'                
                fldr = folder(ns.Experiment& f);
                ff =fullfile(fldr,f.filename);
                if exist(ff,"file")
                    fCntr=fCntr+1;
                    files{fCntr} = ff;
                    info = dir(ff);
                    bytes(fCntr) = info.bytes;                
                else
                    fprintf('File not found: %s\n',ff);
                end
            end
            % Pick the requested ones
            switch sizeOption
                case -1
                    [~,ix] = sort(bytes,'ascend');
                    ix= ix(1);
                case +1
                    [~,ix] = sort(bytes,'descend');
                    ix= ix(1);
                case 0
                    ix = 1:fCntr;
            end
            v =ff(ix);
            if numel(v)==1
                v=v{1};
            end
        end
    end

    methods (Access=protected)
        function makeTuples(tbl,key)
            fldr = folder(ns.Experiment& key);
            ff =fullfile(fldr,key.filename);
            if exist(ff,'file')
                try
                    mv= VideoReader(ff); 
                catch me
                    fprintf('Could not open %s.\n Error: %s\n', ff,me.message)
                    return;
                end
            else
                error('File not found %s',ff);
            end
            tpl = mergestruct(key,struct('nrframes',mv.NumFrames,'width',mv.Width,'height',mv.Height,'framerate',mv.FrameRate));
            insert(tbl,tpl)
        end
    end

end