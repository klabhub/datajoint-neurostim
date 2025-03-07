function   [fc,p,err,src,trg] = pearson(parms,channels)
% Determine FC based on pearson correlation using parameter settings 
% defined in parms for the CChannels in channels.
% 
% If parms.align is defined, data will be retrieved for each of the
% channels using the ns.C/align function (with the align struct used as the
% input parameters for that function). Pearson correlation is then
% determined per trial, and averaged across trials. 
% 
% If parms.align is not defined, Pearson correlation is computed on the 
% entire time course in the experiment. 
% 
% Note that channels can have data from multiple ns.Experiments. With 
% 'align' specified, the data from separate experiments are treated as 
% separate trials. When no align is specified, the time courses of the 
% experiments are concatenated.

if isfield(parms,'align')
    % Use the align function to retrieve all
    % signals
    args =namedargs2cell(parms.align);
    [T,channel] = align(channels,args{:});
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
    % Use the raw signal for the entire experiment(s)
    [X,channel] = fetchn(channels ,'signal','channel');
    % Concatenate channels across experiments 
    [G,channel] = findgroups(channel);
    X = splitapply(@(x) {cat(1,x{:})},X,G);
    X = [X{:}];
    nrChannels = numel(channel);
    % Determine pairwise pearson correlation
    [r,p] = corr(X,"Type","Pearson");
    err = nan(size(r)); % Not defined
end

% Extract unique values; below the diagonal
ix=  reshape(1:(nrChannels*nrChannels),nrChannels,nrChannels);
ix = tril(ix,-1);
ix(ix==0) =[];
fc =r(ix);
p = p(ix);
err = err(ix);
[src,trg] =meshgrid(channel,channel);
src = src(ix); % no distinction between src and trg for Pearson
trg = trg(ix);
end