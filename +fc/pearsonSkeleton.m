function  [srcChannels,trg]  = pearsonSkeleton(parms,srcChannels,trgChannels)
% Function to compute a FC skeleton based on Pearson correlations. 
% This uses the parms.skeleton.upper and .lower values to include only 
% pairs with an FC that falls above the upper percentile or below the
% lower.
%
arguments
    parms (1,1) struct % The parameters to use for FC computation 
    srcChannels (1,1) ns.CChannel % The CChannels to use as sources
    trgChannels (1,1) ns.Channel  % The CChannels to use as targets
end
% Determine pearson FC for all
[FC,p,err,srcChannels,trg] = fc.pearson(parms,srcChannels,trgChannels); %#ok<ASGLU>

% Now determine the skeleton
stay   = false(size(FC)); % By default none
% Include upper percentile
if isfield(parms.skeleton,'upper')
    stay = stay | FC > prctile(FC,parms.skeleton.upper);
end
% Include lower percentile
if isfield(parms.skeleton,'lower')
    stay = stay |  FC < prctile(FC,parms.skeleton.lower);
end
% Return only pairs that the criteria
srcChannels =srcChannels(stay);
trg = trg(stay);