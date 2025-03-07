%{
# Functional Connectivity 
-> ns.C       # Continuous data used to compute FC.
-> fc.Parm    # Parameters that define the FC computation 
---
-> fc.Skeleton = NULL
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

            % Fetch the parms from the FcParms table
            parms = fetch1(ns.FcParm&key,'parms');
            if isa(parms.method,'function_handle')
                % User defined function that takes the key
                % and parms as input and returns a tpl with src, trg,
                % fc, p, and err. Call it
                pairTpl = parms.method(key,parms);
            else
                %  A built-in method.
                switch upper(parms.method)
                    case 'PEARSON'
                        % Pearson correlation between two signals
                        if isfield(parms,'align')
                            % Use the align function to retrieve all
                            % signals
                            args =namedargs2cell(parms.align);
                            [T,channel] = align(ns.C & key,args{:});
                            X= timetableToDouble(T);
                            X =permute(X,[1 3 2]);
                            [~,nrChannels,nrTrials] =size(X);
                            r = zeros(nrChannels,nrChannels,nrTrials);
                            p = zeros(nrChannels,nrChannels,nrTrials);
                            for tr=1:nrTrials
                                [r(:,:,tr),p(:,:,tr)] = corr(X(:,:,tr),"Type","Pearson");
                            end
                            err = std(r,0,3,"omitmissing");
                            r = mean(r,3,"omitmissing");
                            p = 1-mean(p<0.05,3,"omitmissing");
                        else
                            % Use the raw signal for the entire experiment
                            [X,channel] = fetchn(ns.CChannel & key,'signal','channel');
                            nrChannels = numel(channel);
                            % Determine pairwise pearson correlation
                            [r,p] = corr(X,"Type","Pearson");
                            err = nan(size(r));
                        end

                        % Extract below the diagonal
                        ix=  reshape(1:(nrChannels*nrChannels),nrChannels,nrChannels);
                        ix = tril(ix,-1);
                        ix(ix==0) =[];
                        r =r(ix);
                        p = p(ix);
                        err = err(ix);
                        % identify the corresponding src and trg channels.
                        % (not directional in this case, so src could
                        % be labeled trg and vice versa).
                        [srcIx,trgIx] = ind2sub([nrChannels,nrChannels],ix);
                        src =channel(srcIx);
                        trg = channel(trgIx);
                        % Create tpls for each pair
                        pairTpl = struct('fc',num2cell(r(:)), ...
                                         'p',num2cell(p(:)), ...
                                         'err',num2cell(err(:)), ...
                                         'source',num2cell(src(:)), ...
                                         'target',num2cell(trg(:)));
                    otherwise
                        error('Unknown FC methods %s',parms.method)
                end
            end
            % Insert placeholder into the Fc table
            insert(tbl,key);

            % Combine the computed pairwise values with the key
            pairTpl = mergestruct(key,pairTpl);
            % Insert into the FcPair table
            chunkedInsert(ns.FcPair,pairTpl);
        end
    end
end
