%{
# Transformed Epoch. 
-> ns.Epoch
-> ns.TepochParm
---
x : longblob
dependent : varchar(32)
independent : varchar(32) 
%}
classdef Tepoch < dj.Computed & dj.DJInstance   
   methods
       
   end
    methods (Access = protected)
        function makeTuples(tbl, key)
            % Apply a computation/transform to a collection of Epochs and
            % store as Tepoch.

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

            % Compute - uses the ns.cache/compute function  (abstract
            % superclass).
            [T,dv,idv] = compute(ecTbl,parms.fun,parms.options,grouping=parms.grouping,timeWindow=parms.window);
                      
           % Insert in the table
           tpl = key;
           tpl.x = T{1,idv}; % Store the first row only-all others are the same (see fill)
           tpl.dependent = dv;
           tpl.independent = idv;
           insert(tbl,tpl);
           y = T.(dv);
           if isnumeric(y)
               y = num2cell(y,2);
           end
            tpl = mergestruct(key,...
                        struct('y',y,...
                                'channel',num2cell(T.channel),...
                                'trial', num2cell(T.trial)));
           insert(ns.TepochChannel,tpl)

        end
    end

end