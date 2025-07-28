classdef C_layout < dj.DJProperty

    properties (Access = protected)
               
        parent

    end

    properties (Access = protected, Constant)

        identity_vars = {'filename'}
    
    end

    methods

        function self = C_layout(cTbl)
            
            if nargin == 0
                return;
            end
            
            self.parent = cTbl;
            self.make();
            
        end

    end

    methods (Access = protected)                

        function l = get_method(~, key)

            l = getfield(fetch(key, 'info'), 'info');
            if isfield(l, 'layout')
                l = l.layout;
            else
                l = [];
            end

        end

    end

    
end