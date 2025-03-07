%{
# Functional Connectivity skeleton 
-> ns.Session
-> fc.Parm  # Parameters determining how to estimate the skeleton
---
pairs : int     # Number of pairs in the skeleton
%}
% The fc.Parm parms struct should have a member called skeleton with
% .paradigm = the experiments to use for the skeleton
% .method  = The name of a function that computes the skeleton. 
% This function will be called with two inputs.
% The first is the fc.Parm parms structure - it stores the instructions
% on how to compute the skeleton. 
% The second is the ns.CChannel table that should be used to compute the 
% skeleton (it has all channels selected by the combination of .paradigm and 
% session).
% The function should determine the channels based on the .parms and
% ns.CChannel input, and returns these values to the caller.
%
% EXAMPLE:
% fc.peasonSkeleton
% 
% See Also : fc.pearsonSkeleton fc.Fc

classdef Skeleton <dj.Computed
      methods
        function G = summary(tbl)
            T= fetchtable(tbl);
            G = groupsummary(T,["subject" "session_date" "fctag"]);
        end       
    end

    methods (Access=protected)
        function makeTuples(tbl,key)
            parms =fetch1(fc.Parm & key,'parms');
            if isfield(parms,'skeleton')
                % Skeleton instructions provided
                %% Determine which experiments and ns.C rows should be used
                expt = ns.Experiment & key &  in("paradigm",parms.skeleton.paradigm);
                C =     (ns.C & expt & key);              
                assert(exists(ns.CChannel &C),sprintf("No relevant channels available for %s-%s with paradigms %s",key.subject,key.session_date,parms.skeleton.paradigm));

                %% Compute the skeleton
                if ~isempty(which(parms.skeleton.method))
                    % Call the specified function - it should
                    % return source and  target filled with channel numbers.                    
                    channels= feval(parms.skeleton.method,parms,ns.CChannel &C);
                else
                    error("Unknown skeleton method %s",parms.skeleton.method);
                end
                pairTpl = mergestruct(key,struct('channel',num2cell(channels(:))));
                key.pairs = numel(channels);
                insert(tbl,key);
                chunkedInsert(fc.SkeletonChannel,pairTpl);                           
            end
        end
    end
end