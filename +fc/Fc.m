%{
# Functional Connectivity 
-> ns.C         # Continuous data used to compute FC.
-> fc.Skeleton  # Parameters that define the FC computation 
---
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
                tbl (1,1) {mustHaveRows(tbl)}               
                pv.bothSides (1,1) logical % Show symmetric matrix
                pv.alpha (1,1) double = 0.05;                 
            end
            tiledlayout('flow')
            for tpl = fetch(tbl)'
                nexttile
                FC =fcMatrix(tbl&tpl,type =["fc" "p"]);
                FC.fc(FC.p>pv.alpha) = NaN;
                imagesc(FC.fc);
                xlabel 'Channel'
                ylabel 'Channel'                
                title(sprintf('%s for %s-%s@%s.',tpl.fctag,tpl.subject,tpl.session_date,tpl.starttime))                     
            end
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
                [v.channel] =unique([[FC.source];[FC.target]]);
                [~,srcIx] = ismember([FC.source],v.channel);
                [~,trgIx] = ismember([FC.target],v.channel);                
                nrChannel = numel(v.channel);                
                sz = [nrChannel nrChannel];
                for tp = pv.type
                    v.(tp) = nan(sz);
                    v.(tp)(sub2ind(sz,srcIx,trgIx)) = [FC.(tp)];                        
                end
        end
    end

    methods (Access=protected)

        function makeTuples(tbl,key)
            % When calling (par)populate, this function is called with the
            % key containing the information on a set of continuous data
            % (i.e. limited to an experiment), and a skeleton (which is 
            % computed per session, so it can be based on multiple
            % experiments). 
            % In principle FC can be computed for any C data that has
            % multiple channels.

            %% Fetch the parms from the FcParms table
            parms = fetch1(fc.Parm&key,'parms');
            assert(~isempty(which(parms.method)),sprintf('Unknown FC method %s',parms.method))
           
           %% Call the user specified method function 
            % The method takes the parms and channels as input and returns 
            % fc, p, and err, for each src and trg. 
            channels = ns.CChannel & fc.SkeletonChannel & key ; % Restrict channels to the skeleton                        
            [FC,p,err,src,trg] = feval(parms.method,parms,channels);      

            % Insert the key into the Fc table
            % TODO: key.quality = add some generic quality measure...            
            insert(tbl,key);

            % Prepare tuples for the FcPair table
            tpl = struct('fc',num2cell(FC(:)),'p',num2cell(p(:)),'err',num2cell(err(:)),'source',num2cell(src(:)),'target',num2cell(trg(:)));
            tpl = mergestruct(ns.stripToPrimary(fc.Fc,key),tpl);
            % Insert.
            chunkedInsert(fc.FcPair,tpl);
        end
    end
end
