%{
# Functional Connectivity skeleton 
-> ns.Session
-> fc.Parm  # Parameters determining how to estimate the skeleton
-> ns.CChannel
source = 1 : tinyint 
target = 1 : tinyint
---
%}
% The fc.Parm parms struct should have a member called skeleton with 
% .ctag = the C data to use for the skeleton
% .paradigm = the experiments to use for the skeleton
% .method  = how to compute the skeleton. 
% The .method can be ALL2ALL to include all channels in the skeleton, or 
% it can be the name of a function that computes the skeleton. This
% function will be called with three inputs. The first is the Skeleton tuple
% containing the primary keys and empty source and target. 
% The second is the fc.Parm.parms structure - it should store the instructions 
% on how to compute the skeleton. The third is the ns.CChannel table that 
% should be used to compute the skeleton. 
% The fun computes the source and target channels based on the .parms and 
% ns.CChannel input, sets the .source  and .target fields of the Skeleton tuple
% and returns this to the caller. 
% 
% See Also : fc.pearsonSkeleton fc.Fc
               
classdef Skeleton <dj.Computed
    properties (Dependent)
        keySource
        
    end

    methods
        function G = summary(tbl)
                T= fetchtable(tbl);
                G = groupsummary(T,["subject" "session_date" "fctag"]);
        end
        function v = get.keySource(~)
            %  key source is just the session and parm
            v = ns.Session *fc.Parm;
        end
    end

    methods (Access=protected)
        function makeTuples(tbl,key)
            parms =fetch1(fc.Parm & key,'parms');
            if isfield(parms,'skeleton')
                % Skeleton instructions provided
                %% Determine which experiments and ns.C rows should be used
                expt = ns.Experiment & key &  in("paradigm",parms.skeleton.paradigm);
                C =     (ns.C & expt) & sprintf('ctag="%s"',parms.skeleton.ctag);                
                allChannels = fetch(ns.CChannel & C ,'channel');
                %% Compute the skeleton
                switch upper(parms.skeleton.method)
                    case "ALL2ALL"
                        src =  [allChannels.channel];
                        trg = [allChannels.channel]; % Same
                    otherwise
                        if ~isempty(which(parms.skeleton.method))
                            % Call the specified function - it should
                            % return source and  target filled with channel numbers.
                            [src,trg]= feval(parms.skeleton.method,parms,ns.CChannel &C);
                        else
                            error("Unknown skeleton method %s",parms.skeleton.method);
                        end
                end
                
                %% Determine source and target - usually channels are both
                bothSrcAndTrg = intersect(src,trg);               
                tpl = allChannels(ismember([allChannels.channel],bothSrcAndTrg));
                [tpl.source] =deal(1);
                [tpl.target] =deal(1);                    
                onlySrc= setdiff(src,trg);
                onlySrcTpl = allChannels(ismember([allChannels.channel],onlySrc)); 
                if ~isempty(onlySrcTpl)
                    [onlySrcTpl.source] = deal(1);
                    [onlySrcTpl.target] = deal(0);
                    tpl = [tpl;onlySrcTpl];
                end
                onlyTrg= setdiff(trg,src);
                onlyTrgTpl = allChannels(ismember([allChannels.channel],onlyTrg)); 
                if ~isempty(onlyTrgTpl)
                    [onlyTrgTpl.source] = deal(0);
                    [onlyTrgTpl.target] = deal(1);
                    tpl = [tpl;onlyTrgTpl];
                end
                
                [tpl.fctag] = deal(key.fctag);
                insert(tbl,tpl);
            %else- no skeleton instructions - no entry created. FC will be
            %computed on all channels.
            end
        end
    end
end