%{
# TepochParm: parameters to transform epoch data.
ttag : varchar(32)  # tag for this transformed epoch
etag : varchar(32)  # which epochs to transform
---
description = NULL : varchar(512) # Brief description.
fun : longblob       # Struct that will be passed as fun to ns.cache.compute
parms : longblob     # Struct of other name-value inputs for ns.cache.compute
%}
%
classdef TepochParm < dj.Lookup & dj.DJInstance
    methods
        function insert(tbl,tpl)
            arguments
                tbl (1,1)
                tpl (:,1) struct
            end
            tuples = tpl;
            for i = 1:numel(tuples)
                pv = namedargs2cell(tuples(i));
                tuple = set_defaults(pv{:});
                if i == 1
                    tpl = repmat(tuple,size(tuples));
                end
                tpl(i) = tuple;
            end
            tpl = makeMymSafe(tpl);
            insert@dj.Lookup(tbl,tpl);
        end
    end
end

function tpl = set_defaults(tpl)
% Mainly serves to check the correct data types.
arguments
    tpl.ttag (1,1) string 
    tpl.etag (1,1) string 
    tpl.description (1,1) string = ""
    tpl.parms (1,1) struct = struct('channel',[],'trial',[],'timeWindow',[-inf inf],'average',string.empty);
    tpl.fun (1,1) struct
end
            

end