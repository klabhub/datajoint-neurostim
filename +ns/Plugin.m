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
            % This is called by updateWitFileContents from ns.Experiment\

            key.plugin_name = plg.name;
            insert(self,key);
            make(ns.PluginParameter,key,plg.prms);

            % CiC has some properties that are useful to store, but they
            % are not regular properties (and therefore not added
            % automatically via the pluginparameter make function). We
            % handle those separately here.
            if strcmpi(plg.name,'cic')
                % 1.  create a blockName property that stores the names of
                % the blocks.
                blockNameTpl = fetch(ns.PluginParameter & key & 'property_name="block"','*');
                blockNameTpl.property_name = 'blockName';
                blockNr  = blockNameTpl.property_value ;
                blockNameTpl.property_value = cell(size(blockNr));
                stay =blockNr>0;
                blockNames = {plg.blocks(blockNr(stay)).name};
                [blockNameTpl.property_value{stay}] =deal(blockNames{:});
                insert(ns.PluginParameter,blockNameTpl);
                % 2. More can be added here
            end
        end

        function what(tbl)
            % Function to show the parameters for each plugin in the table.
            % Includes links to retrieve the parameter values.
            expt= fetch(tbl,'subject','session_date','starttime','LIMIT 1');
            fprintf("PluginParameter values are taken from %s:%s:%s\n",expt.subject,expt.session_date,expt.starttime)
            for plg = fetch(tbl)'
                fprintf(2,'\t%s:\n',plg.plugin_name)
                for tp=["Global" "Parameter" "Event" "Bytestream"]
                    T = fetchtable(ns.PluginParameter & plg & struct('property_type',tp));
                    if ~isempty(T)
                        line ='';
                        if tp=="Parameter"
                            extra = ",'AtTrialTime',0";
                        else
                            extra = "";
                        end
                        for prm = 1:height(T)
                            cmd = sprintf('fprintf(2,''%s'');%s = get(ns.Experiment & struct(''subject'',''%s'',''session_date'',''%s'',''starttime'',''%s''),''%s'',''prm'',''%s'',''what'',''data''%s)',plg.plugin_name,T.property_name(prm),T.subject(prm),T.session_date(prm),T.starttime(prm),T.plugin_name(prm),T.property_name(prm),extra);
                            line= [line ' ' sprintf('<a href="matlab:%s">%s</a>',cmd,T.property_name(prm))]; %#ok<AGROW>
                        end
                        disp(sprintf('\t\t%-10s: %s\t',tp, line)) %#ok<DSPSP>
                    end
                end
            end
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

