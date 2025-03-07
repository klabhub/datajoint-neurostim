%{
# Functional Connectivity 
-> ns.C       # Continuous data used to compute FC.
-> fc.Parm    # Parameters that define the FC computation 
---
pairs           : int  # Number of pairs in the FC
quality = NULL : float # Some measure of quality of the FC (optional)
%}
%
% Functional connectivity table. 
%  The ns.C foreign key identifies the continuous data for which functional
%  connectivity is computed.
% The fc.Parm is a lookup table that defines how FC is computed. The user
% has full control by specifying a method in the fc.Parm parms. 
%
% The dj.Part table fc.FcPair stores the actual FC values per pair
%
% See Also fc.Parm
classdef Fc < dj.Computed
     methods
        function plot(tbl,pv)
            % Rudimentary plot function. Needs to be extended with options
            arguments
                tbl (1,1) {mustHaveRows(tbl,1)}               
                pv.bothSides (1,1) logical % Show symmetric matrix
                pv.alpha (1,1) double = 0.05;                 
            end
            FC =fcMatrix(tbl,type =["fc" "p"]);
            FC.fc (FC.p<pv.alpha) = NaN;
            imagesc(FC.fc);
            xlabel 'Target'
            ylabel 'Src'
            key =fetch(tbl);
            title(sprintf('%s for %s-%s@%s.',key.fctag,key.subject,key.session_date,key.starttime))                     
        end

        function v = fcMatrix(tbl,pv)
            % Convenience function that returns a struct with FC as a
            % matrix, together with the channels that the FC was
            % calculated for and, optionally, the associated p and err
            % values form the FcPair table.
            arguments
                tbl (1,1) fc.Fc
                pv.type (1,:) string = "fc"  % Set to ["fc" "p" "err" ]  to include p and err values.
            end
                columns =cellstr([pv.type "source" "target"]);
                FC = fetch(fc.FcPair & tbl,columns{:});
                [v.src,~,srcIx] =unique([FC.source]);
                [v.trg,~,trgIx] =unique([FC.target]);
                nrSource = numel(v.src);
                nrTarget = numel(v.trg);
                for tp = pv.type
                    m = nan(nrSource,nrTarget);
                    m(sub2ind([nrSource nrTarget],srcIx,trgIx)) = [FC.(tp)];    
                    v.(tp) = m;
                end
        end
    end

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
            assert(exists(skeleton),"No skeleton found for %s-%s@%s. Populate fc.Skeleton first?.",key.subject,key.session_date,key.starttime)                
            % Pass the CChannel to the method, together with issource/istarget
            % information from the skeleton (in case the method deals with
            % these differently).
            
            % The function takes the parms and channels as input and returns 
            % fc, p, and err, for each src and trg. 
            [FC,p,err,source,target] = feval(parms.method,parms,ns.CChannel * proj(skeleton,'issource','istarget'));
            key.pairs = numel(source);
            % Insert the result into the Fc table
            insert(tbl,key);

            % Prepare tuples for the Pair table
            tpl = struct('fc',num2cell(FC(:)),'p',num2cell(p(:)),'err',num2cell(err(:)),'source',num2cell(source(:)),'target',num2cell(target(:)));
            tpl = mergestruct(ns.stripToPrimary(fc.Fc,key),tpl);
            chunkedInsert(fc.FcPair,tpl);
        end
    end
end
