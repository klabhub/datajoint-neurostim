classdef VariableLengthCache < ns.cache
    % Database-free cache fixture for heterogeneous result widths.
    methods
        function self = VariableLengthCache(T)
            self.T = T;
            self.samplingRate = 1;
        end
    end
    methods (Access=protected)
        function fill(~)
            % The constructor supplies the already-filled cache table.
        end
        function src = getCacheQuery(~)
            src = [];
        end
    end
end