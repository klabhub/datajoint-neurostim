%{
# TepochParm: parameters to transform epoch data.
ttag : varchar(32)  # tag for this transformed epoch
etag : varchar(32) # which epochs to transform 
---
fun : longblob # Function to do the transform (see ns.cache/compute)
window  = NULL  : tinyblob          # Start and stop time of the epoch. Defaults to entire epoch
channels = NULL  : blob              # Channels to include. Defaults to all in the etag
average = NULL :blob               # Things to avere over  one ore more of these ["trial" "channel"];           
trials  = NULL : blob               # Trials to include. Defaults to all in the etag    
%}


classdef TepochParm < dj.Lookup & dj.DJInstance
    methods 
    function insert(tbl,tpl)
        arguments
            tbl (1,1)
            tpl (:,1) struct
        end

    % TODO define insert to check valid values
%               assert(isempty(parms.average) || all(ismember(parms.average,["trial" "channel"])),"Tepoch averaging must ...")

        tpl = makeMymSafe(tpl);
        insert@dj.Lookup(tbl,tpl);
    end
    end

end