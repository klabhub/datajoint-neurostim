%{
# Functional Connectivity 
-> ns.C       # Continuous data used to compute FC.
-> fc.Parm    # Parameters that define the FC computation 
---
-> fc.Skeleton 
%}
%
% Functional connectivity table. Draft version.
%  The ns.C foreign key identifies the continuous data for which functional
%  connectivity is computed.
% The fc.Oarm is a lookup table that defines how FC is computed. The user
% has full control by specifying a function handle as a method in the
% fc.Parm parms. 
% 
% Currently this table only serves as bookkeeping (that FC has been
% computed); it could be extended with additional properties that apply to
% the entire FC (properties that are not specific to a pair, and not
% already in the fc.Parm). 
%
% The part table fc.Pair stores the actual values plus stats of the computed 
% connectivity between pairs. This too could be extended with additional properties. 
% 
classdef Fc < dj.Computed
    methods (Access=protected)

        function makeTuples(tbl,key)
            % When calling (par)populate, this function is called with the
            % key containing the information on a set of continuous data.
            % In principle FC can be computed for any C data that has
            % multiple channels.

            %% Fetch the parms from the FcParms table
            parms = fetch1(fc.Parm&key,'parms');
            assert(~isempty(which(parms.method)),sprintf('Unknown FC method %s',parms.method))
            
            %% Determine which channels are in the skeleton
            skeleton =  fc.Skeleton & key;

            % the function takes the parms and channels as input and returns 
            % fc, p, and err, for each src and trg. 
            [FC,p,err,src,trg] = feval(parms.method,parms,ns.CChannel & skeleton);
            
                
                % Insert placeholder into the Fc table
            insert(tbl,key);

            % Combine the computed pairwise values with the key
            pairTpl = mergestruct(key,pairTpl);
            % Insert into the FcPair table
            chunkedInsert(ns.FcPair,pairTpl);
        end
    end
end
