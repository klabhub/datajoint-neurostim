%{
# List of paradigms for which a leadfield should be computed.
->ns.LeadfieldParm 
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
                pv.description (1,1) string
                pv.parms  (1,1) struct 
            end
            
             % The mode can be 
             %  e2d (loop over electrodes)
             %  d2e (loop over dipoles), 
             %  auto (guess which loop is faste and use that), 
             %  ft to use the Fieldtrip method (equivalent to e2d)
             % Note that the Fieldtrip method is serial, while the other
             % modes can run in a parfor (see nsParPool how to setup the
             % pool, if no pool is availble, they run serially. The
             % resulting leadfield is the same for all modes only the
             % computation time varies.
            defaults.mode = 'auto' ;
            defaults.cfg.reducerank = 'no'; % The parameters that are passed to ft_prepare_leadfield
            defaults.cfg.backproject = 'yes';
            defaults.cfg.normalize = 'no';
            defaults.cfg.normalizeparam = 0.5;
            defaults.cfg.weight = [];

            pv.parms= mergeDefaults(pv.parms,defaults);
            
           
        end
    end



end