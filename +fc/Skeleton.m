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
% fc.pearsonSkeleton
% 
% See Also : fc.pearsonSkeleton, fc.Fc

classdef Skeleton < dj.Computed
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
                assert(exists(ns.CChannel & C),sprintf("No relevant channels available for %s-%s with paradigms %s",key.subject,key.session_date,parms.skeleton.paradigm));

                %% Compute the skeleton
                if ~isempty(which(parms.skeleton.method))
                    % Call the specified function - it should
                    % return source and  target filled with channel numbers.
                    channels = ns.CChannel & C;
                    
                    % define the kind of dataset that will be created
                    if isfield(parms,'align')
                        % bad_rois is informative only,
                        % maybe useful for preprocessing checks. Where to
                        % save?
                        [D,channelInMatrix,bad_rois] = fc.alignForFC(parms,channels); 
                    end
                    % compute FC method selected in parms.skeleton.method
                    [FC,p,err,src,trg] = feval(parms.skeleton.method,parms,D,channelInMatrix); %#ok<ASGLU>
                    % % Average per channel
                    [G,allChannels] = findgroups([src;trg]);
                    % nrChannels = numel(allChannels);
                    % mFC =accumarray(G',FC',[nrChannels 1],@mean);
                    % % Then find the upper and lower percentiles.    
                    % stay = false(nrChannels,1);
                    % % Include upper percentile
                    % if isfield(parms.skeleton,'upper')
                    %     stay = stay | mFC > prctile(mFC,parms.skeleton.upper);
                    % end
                    % % Include lower percentile
                    % if isfield(parms.skeleton,'lower')
                    %     stay = stay |  mFC < prctile(mFC,parms.skeleton.lower);
                    % end
                    % % Return only channel numbers that meet the criteria
                    % channels =allChannels(stay);
                    channels = allChannels;

                else
                    error("Unknown skeleton method %s",parms.skeleton.method);
                end
                % add the table fc.SkeletonChannel
                pairTpl = mergestruct(key,struct('channel',num2cell(channels(:))));
                % TODO: this is incorrectly named. Should be channels or rois.
                key.pairs = numel(channels);
                insert(tbl,key);
                chunkedInsert(fc.SkeletonChannel,pairTpl);

                % add fc.SkeletonFcPair
                tpl = struct('fc',num2cell(FC(:)),'p',num2cell(p(:)),'err',num2cell(err(:)),'source',num2cell(src(:)),'target',num2cell(trg(:)));
                tpl = mergestruct(ns.stripToPrimary(fc.Skeleton,key),tpl);
                chunkedInsert(fc.SkeletonFcPair,tpl);

            end
        end
    end
end
