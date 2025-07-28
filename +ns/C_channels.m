classdef C_channels < dj.DJProperty

    properties (Access = protected)
               
        parent

    end

    properties (Access =protected, Constant)

        identity_vars = {'filename'}
    
    end

    methods

        function self = C_channels(cTbl)

            if nargin == 0
                return
            end
            
            self.parent = cTbl;
            self.make();
            
        end

    end

    methods (Access = protected)                

        function ch = get_method(self, key)

            ch = fetch(ns.CChannel & self.parent & key, 'channel');
            ch = [ch.channel];

        end

    end

    
end