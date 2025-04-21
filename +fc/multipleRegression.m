function [fc,p,err,src,trg] = multipleRegression(parms,D,channelInMatrix,sk)
    % for each target, get all its sources that have non-zero fc.
    % define a regression target ~ B*sources + e, use estimated B as new
    % entries in the FcPair tuples.


    T  = fetchtable(sk,'source','target','fc');
    graph = table2array(T(:,["source","target","fc"]));
    
    nrChannels = numel(channelInMatrix);
    fc = zeros(nrChannels,nrChannels);
    % empty for now, what to save here, if any?
    p =  zeros(nrChannels,nrChannels);
    err = zeros(nrChannels,nrChannels);
    
    % loop for each target
    for r = 1:nrChannels
        trg = channelInMatrix(r);
        % get all adjacencies containing trg, with non-zero fc (in the skeleton)
        idx = (graph(:,1)==trg | graph(:,2)==trg) & graph(:,3) ~= 0;
        pairs =  graph(idx,1:2);
        adj =  pairs(:);
        % remove target channel
        adj(adj == trg) = [];
        % "unique" is just an extra precaution (probably not needed)
        adj = sort(unique(adj));
        % get the data for y and X
        y = D(:,channelInMatrix == trg);
        y = zscore(y);
        [~,src_idx] = ismember(adj,channelInMatrix);
        X = D(:,src_idx);
        X = zscore(X,0,1);
        B = regress(y,X);
        fc(r,src_idx) = B;
    end
    
    % Transform FC data into a list of tuples:
    % source(src),target(trg),fc-value(fc)
    % symmetrize the fc weights by taking the average of (i,j) and (j,i)
    fc = (fc + fc')/2;
    % Extract unique values; below the diagonal
    ix = find(tril(true(nrChannels), -1));
    % include also zero values, we  want all entries
    % to make an FC matrix reconstruction easier to visualize.
    fc = fc(ix);
    % empty for now, not sure what to add here.
    p = p(ix);
    err = err(ix);
    % sources and targets src -- trg
    [src,trg] = meshgrid(channelInMatrix,channelInMatrix);
    src = src(ix);
    trg = trg(ix);

end
