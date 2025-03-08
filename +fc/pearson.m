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
arguments
    parms (1,1) struct
    channels (1,1)  ns.CChannel
end

if isfield(parms,'align')
    % Use the ns.C/align function to retrieve specific
    % signals from the CChannels.
    args =namedargs2cell(parms.align);  % All parms members are passed to align
    allC =ns.C &channels;
    warning('off', 'DataJoint:longCondition')
    fc = [];
    p = [];
    for thisC= fetch(allC)'
        % Loop over rows in C ; they correspond to experiments
        [T,~,channelInMatrix] = align(allC & thisC ,'channel',proj(channels &thisC), args{:});
        X= timetableToDouble(T);
        % Sort the channels to make sure src/trg pairs are all in the upper
        % right corner of a matrix. (For consistency with the method
        % below and ease of visualization in fc.Fc/Plot).
        [channelInMatrix,sorted] = sort(channelInMatrix);
        X = X(:,:,sorted);        
        X =permute(X,[1 3 2]); % Put trials last
        [~,nrChannels,nrTrials] =size(X);
        thisFc = zeros(nrChannels,nrChannels,nrTrials);
        thisP = zeros(nrChannels,nrChannels,nrTrials);
        % FC per trial
        for tr=1:nrTrials
            [thisFc(:,:,tr),thisP(:,:,tr)] = corr(X(:,:,tr),"Type","Pearson");
        end
        fc = cat(3,fc,thisFc); % Concatenate trials across experiments
        p  = cat(3,p,thisP);
    end
    % After collating across trials and experiments (the latter only for
    % skeleton computation), take the average.
    err = std(fc,0,3,"omitmissing");  % Err is the stdev of the FC.
    fc= mean(fc,3,"omitmissing");     % FC is the mean across trials
    p= mean(p<0.05,3,"omitmissing"); % p is the fraction of significant ones
    warning('on', 'DataJoint:longCondition')
else
    % Use the raw signal for the entire experiment(s)
    [X,channelInMatrix] = fetchn(channels,'signal','channel','ORDER BY channel');
    % Concatenate channels (and across experiments for skeletons).
    [G,channelInMatrix] = findgroups(channelInMatrix);
    X = splitapply(@(x) {cat(1,x{:})},X,G);
    X = [X{:}];
    % X now contains the signal for each channel. Channel numbers are
    % defined in channelInMatrix
    [fc,p,rl, ru] = corrcoef(X);
    err = ru-rl;
end

% Extract unique values; below the diagonal
nrChannels = size(fc,1);
ix=  reshape(1:(nrChannels*nrChannels),nrChannels,nrChannels);
ix = tril(ix,-1);
ix(ix==0) =[];
fc=fc(ix);
p= p(ix);
err = err(ix);
[src,trg] =meshgrid(channelInMatrix,channelInMatrix);
src = src(ix);
trg = trg(ix);

end