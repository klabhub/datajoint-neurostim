%{
# A class representing a Neurostim Plugin.
plugin_name : varchar(25)  # The name of the plugin
-> ns.Experiment           # The corresponding experiment (FK)
---
%}
% 
% BK  - April 2022
classdef Plugin < dj.Manual
    methods  (Access=public) % ns.scan needs access
        function key= make(self,key,plg)
            % function key= make(self,key,plg)
            % This is called by updateWitFileContents from ns.Experiment
            
            key.plugin_name = plg.name;
            insert(self,key);
            make(ns.PluginParameter,key,plg.prms);   
           
        end 
    end

    methods (Access=?ns.Experiment)
        function nwbRoot = nwb(tbl,nwbRoot,pv)
            arguments
                tbl (1,1) ns.Plugin 
                nwbRoot (1,1) NwbFile  % The root element
                pv (1,1) struct % THe struct of pv pairs set in ns.Experiment/nwbExport
            end
            %Export to NWB format. Called from ns.Experiment/nwbExport only
            for plgTpl = fetch(tbl,'*')'               
                get(ns.PluginParameter&plgTpl,nwbRoot,plgTpl.plugin_name);                    
            end
        end

    end
end

