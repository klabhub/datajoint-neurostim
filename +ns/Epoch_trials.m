classdef Epoch_trials < dj.DJProperty

    properties (Access = protected)
               
        parent

    end

    properties (Access = protected, Constant)

        identity_vars = {'filename', 'etag'}
    
    end

    methods

        function self = Epoch_trials(eTbl)
            
            if nargin == 0
                return
            end
            
            self.parent = eTbl;
            self.make();
            
        end

    end

    methods (Access = protected)                

        function t = get_method(~, key)

            dTbl = ns.Dimension & key;
            
            t = fetch(key, 'time');
            t = t.time;

            if numel(t)==3
                t = linspace(t(1),t(2),t(3))';
            end

        end

    end

    
end