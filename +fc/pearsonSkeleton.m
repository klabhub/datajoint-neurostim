function  [channelsInSkeleton]  = pearsonSkeleton(parms,channels)
% Function to compute a FC skeleton based on Pearson correlations. 
% This uses the parms.skeleton.upper and .lower values to include only 
% pairs with an average FC that falls above the upper percentile or below the
% lower.
% EXAMPLE
% parms.upper = 90
% parms.lower  =10
% takes the lower and upper 10th percentile (e.g. least connected (which could be negative
% and most positive).
% parms.upper = 100 
% parms.lower = 10  
% would only take the least connected. 

arguments
    parms (1,1) struct % The parameters to use for FC computation 
    channels (1,1) ns.CChannel % The CChannels to use as sources    
end
% Determine pearson FC for all channels
[FC,p,err,src,trg] = fc.pearson(parms,channels); %#ok<ASGLU>
% Average per channel
[G,allChannels] = findgroups(src);
nrChannels = numel(allChannels);
mFC =accumarray(G',FC',[nrChannels 1],@mean);
% Then find the upper and lower percentiles.    
stay = false(nrChannels,1);
% Include upper percentile
if isfield(parms.skeleton,'upper')
    stay = stay | mFC > prctile(mFC,parms.skeleton.upper);
end
% Include lower percentile
if isfield(parms.skeleton,'lower')
    stay = stay |  mFC < prctile(mFC,parms.skeleton.lower);
end
% Return only channel numbers that meet the criteria
channelsInSkeleton =allChannels(stay);
