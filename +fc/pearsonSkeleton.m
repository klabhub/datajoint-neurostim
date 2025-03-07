function  [src,trg]  = pearsonSkeleton(parms,channels)
% Function to compute a FC skeleton based on Pearson correlations. 
% This uses the parms.skeleton.upper and .lower values to include only 
% pairs with an FC that falls above the upper percentile or below the
% lower.
%
arguments
    parms (1,1) struct % The parameters to use for FC computation 
    channels (1,1) ns.CChannel % The CChannels to use as sources    
end
% Determine pearson FC for all channels
[FC,p,err,src,trg] = fc.pearson(parms,channels); %#ok<ASGLU>

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
% Return only src and trg channel numbers that meet the criteria
src =src(stay);
trg = trg(stay);