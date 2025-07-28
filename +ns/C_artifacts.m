classdef C_artifacts < dj.DJProperty

    properties (Access = protected)
               
        parent

    end

    properties (Access =protected, Constant)

        identity_vars = {'filename', 'ctag'}
    
    end

    methods

        function self = C_artifacts(cTbl)

            if nargin == 0
                return
            end
            
            self.parent = cTbl;
            self.make();
            
        end

    end

    methods (Access = protected)                

        function art = get_method(~, key)

            info = fetch(key, 'info');
            if isfield(info.info, 'noisyChannels')

                art = info.info.noisyChannels;

            else

                art = struct.empty();

            end

        end

    end

    
end