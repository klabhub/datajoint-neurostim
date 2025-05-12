classdef C_channel_coordinates < dj.DJProperty

    properties (Access = protected)
               
        parent

    end

    properties (Access = protected, Constant)

        identity_vars = {'filename'}
    
    end

    methods

        function self = C_channel_coordinates(cTbl)
            
            if nargin == 0
                return;
            end
            
            self.parent = cTbl;
            self.make();
            
        end

    end

    methods (Access = protected)                

        function coord = get_method(~, key)

            ccTbl = fetch(ns.CChannel & key,'channelinfo');
            info = [ccTbl.channelinfo];
            coord = [info.X; info.Y; info.Z]';

        end

    end

    
end