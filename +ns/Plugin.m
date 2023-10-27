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
end

