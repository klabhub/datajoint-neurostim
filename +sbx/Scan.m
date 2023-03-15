%{
# A scan associated with this experiment
->ns.Experiment
---
info : longblob # The complete info struct
nrframes : int # Number of frames in the scan
%}

classdef Scan < dj.Computed
    properties (Dependent)
        keySource
    end
    methods
        function v = get.keySource(~)
            % Restrict to experiments that have an sbx file
            v = ns.Experiment & (ns.File & 'extension=''.sbx''');
        end
    end
    methods (Access =protected)
        function makeTuples(tbl,key)
            % Look for .sbx file in the ns.File table
            v = ns.File & key & 'extension=''.sbx''';
            if exists(v)
                fname = fetch1(v,'filename');
                fldr = folder(ns.Experiment & key);
                ff =fldr + fname(1:end-3) + 'mat';
                if exist(ff,'file')
                    load(ff ,'info'); % Get the info struct
                else
                    error('File %s not found',strrep(ff,'\','/'))
                end
                N= sbx.nrFrames(fldr + fname,info);
                % Insert in the table
                tpl = mergestruct(key,struct('info',info,'nrframes',N));
                insert(tbl,tpl);
            end
        end
    end
end