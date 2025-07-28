classdef C_time < dj.DJProperty

    properties (Access = protected)
               
        parent

    end

    properties (Access = protected, Constant)

        identity_vars = {'filename'}
    
    end

    methods

        function self = C_time(cTbl)
            
            if nargin == 0
                return
            end
            
            self.parent = cTbl;
            self.make();
            
        end

    end

    methods (Access = protected)                

        function t = get_method(~, key)

            t = fetch(key, 'time');
            t = t.time;

            if numel(t)==3
                t = linspace(t(1),t(2),t(3))';
            end

        end

    end

    
end