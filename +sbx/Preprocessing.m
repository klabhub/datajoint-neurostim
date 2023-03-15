%{
#  Preprocessing parameter sets used for calcium imaging data
    name:  varchar(128)
    ---
    toolbox : enum('suite2p','caiman')
    description : varchar(1024)
    parms : longblob  # struct of all parameters
%}
classdef Preprocessing < dj.Lookup


    methods (Static)
        function default(tbx)
            % Create a preprocessing option called 'default' that simply uses the toolbox defaults 
            switch (tbx)
                case 'suite2p'                        
                        tpl = struct('parms',struct,'name','default','toolbox','suite2p','description','Toolbox defaults');                      
                        insert(sbx.Preprocessing,tpl);
                case 'caiman'
                %otherwise
                % Not possible
            end
        end
    end
end