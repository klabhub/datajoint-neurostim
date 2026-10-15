classdef CacheInsertValidator < ns.cache
    % Exercise the shared insertion guard without a database connection.
    properties (Constant)
        header = struct('names',{{'signal'}})
    end
    methods (Access=protected)
        function src = getCacheQuery(~)
            src = [];
        end
    end
end
