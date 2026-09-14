%{
# TepochParm: parameters to transform epoch data.
ttag : varchar(32)  # tag for this transformed epoch
etag : varchar(32)  # which epochs to transform
---
fun : longblob       # Struct that will be passed as fun to ns.cache.comput
parms : longblob     # Struct of other name-value inputs for ns.cache.compute
%}

classdef TepochParm < dj.Lookup & dj.DJInstance
    methods
        function insert(tbl,tpl)
            arguments
                tbl (1,1)
                tpl (:,1) struct
            end

            tpl = namedargs2cell(tpl);
            tpl = set_defaults(tpl{:});
            tpl = makeMymSafe(tpl);
            insert@dj.Lookup(tbl,tpl);
        end
    end
end

function tpl = set_defaults(tpl)
arguments
    tpl.ttag char
    tpl.etag char
    tpl.parms (1,:) struct = struct('channel',[],'trial',[],'timewindow',[-inf inf],'average',"");
    tpl.fun (1,1) struct
end
            

end