%{
# Transformed Epoch. 
-> ns.Epoch
-> ns.TepochParm
---
x : longblob
%}
classdef Tepoch < dj.Computed & dj.DJInstance   
   methods
       
   end
    methods (Access = protected)
        function makeTuples(tbl, key)
            parms = fetch(ns.TepochParm & key,'*');
            

            % Restrict the epoch channels and trials if requested
            ecTbl = ns.EpochChannel & key;
            if ~isempty(parms.channels)
                 ecTbl = ecTbl & struct('channel',num2cell(parms.channels)); 
            end
            if ~isempty(parms.trial)
                 ecTbl = ecTbl & struct('trial',num2cell(parms.trial)); 
            end            
            % Group if requested
            if isempty(parms.grouping)
                parms.grouping = ["subject" "session_date" "starttime" "condition" "trial" "channel"];
            end
            % Restrict to a specific window
            if isempty(parms.window)                
                parms.window = fetch(ns.EpochParm & key,'window');                
            end

            ecTbl.cache; % Fill the cache with the latest
            % Compute 
            [T,dv,iv] = compute(ecTbl,parms.fun,parms.options,grouping=parms.grouping,timeWindow=parms.window);
                      
           % Insert in the table
           tpl = key;
           tpl.x = T.(iv);
           insert(tbl,tpl);
            tpl = mergestruct(key,...
                        struct('y',T.(dv),...
                                'channel',num2cell(T.channel),...
                                'trial', num2cell(T.trial)));
           insert(ns.TepochChannel,tpl)

        end
    end

end