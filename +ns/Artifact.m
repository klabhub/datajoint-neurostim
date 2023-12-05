%{
# Identifies trials and time periods that should not be used for analysis
-> ns.C
-> ns.ArtifactParm
--- 
trial =NULL :      longblob   # A vector of trials that are considered artifacts and should be removed from analysis
start = NULL : longblob  # A vector of neurostim times that indicate the start of an artifact period
stop =NULL  : longblob   # A vector of neurostim times that indicate the end of an artifact period
%
% This table is used by ns.C/align to remove trials or time periods. To
% setup artifact detection, add an entry to ns.ArtifactParm that defines
% the function that identifies the trials/timepoints that contain artifacts.
% 
%}

classdef Artifact < dj.Computed
    properties (Dependent)
        keySource
    end

    methods
        function v = get.keySource(~)
            % Restricted to ns.C tuples with the ctag specified in
            % DirtyParm 
            allTpl = [];                       
            for thisPrm= fetch(ns.ArtifactParm,'c')'
                tbl = (ns.C*ns.CParm) & struct('ctag',thisPrm.c);    
                % Would like to concatenate this tbl with the next row but
                % this does not work with the | operator. Instead, concatenate
                % tuples of primary keys
                thisTpl = fetch(tbl);
                if isempty(allTpl)
                    allTpl = thisTpl;
                else
                    allTpl  = catstruct(1,allTpl,thisTpl);
                end
            end
            % And then restrict the full table by the set of found tuples.
            v = (ns.C*proj(ns.ArtifactParm,'atag','fun','description','parms')) & allTpl;
        end 
    end
    methods (Access=protected)
        function makeTuples(tbl,key)
            prms = fetch(ns.ArtifactParm &key,'*');
            assert(~isempty(which(prms.fun)),'Artifact detection function %s not found.',prms.fun);
            [experiment,channel] = feval(prms.fun,ns.C&key,prms.parms); % Pass to handler
            % Setup the artifact table (applies to all channels)
            aKey = mergestruct(key,experiment);
            insert(tbl,aKey);
            if ~isempty(channel)
                % Add channel specific artifact info
                acKey = mergestruct(key,channel);
                insert(ns.ArtifactChannel,acKey);
            end
        end
     end

end