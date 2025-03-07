%{ 
# Functional Connectivity skeleton 
-> ns.Session
-> fc.Parm  # Parameters determining how to estimate the skeleton
---
source : blob                 #
target: blob
%}

classdef Skeleton <dj.Computed

    methods (Access=protected)
        function makeTuples(tbl,key)
            parms =fetch1(ns.FcParm & key,'parms');
            if isfield(parms,'skeleton')
                % Skeleton instructions provided 
                % This should have
                % .ctag = the C data to use for the skeleton
                % .paradigm = the experiments to use for the skeleton

                % Determine which ns.C rows should be used

                ns.C & key & "ctag=" + parms.ctag & (ns.Experiment & )

        end
    end
end