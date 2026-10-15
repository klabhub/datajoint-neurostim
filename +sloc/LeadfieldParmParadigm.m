%{
# List of paradigms for which a leadfield should be computed.
->sloc.LeadfieldParm
->ns.Paradigm
%}

classdef LeadfieldParmParadigm < dj.Lookup & dj.DJInstance
     methods
        function insert(self, tuples, varargin)
           % Overload the insert method to do argument validation and setup
           % defaults before insertion
            pv = namedargs2cell(tuples);
            tuples = self.validate(pv{:});
            % Call the superclass function after validation
            insert@dj.Lookup(self,makeMymSafe( tuples), varargin{:})
        end
    end

    methods (Static, Access = protected)
        function pv = validate(pv)
            arguments
                pv.lftag (1,1) string
                pv.name (1,1) string
            end
            
           
        end
    end



end