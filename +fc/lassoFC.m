function [fc,p,err,src,trg] = lassoFC(parms,D,channelInMatrix)


    % compute lasso FC with D
    n_rois = size(D,2);
    fc = zeros(n_rois,n_rois);
    % empty for now, what to save here, if any?
    % maybe NaN to initialize
    p =  zeros(n_rois,n_rois);
    err = zeros(n_rois,n_rois);
    warning('off' , 'stats:lasso:StandardizeIgnored');
    % if lambda is cv, select lambda as the average of 20 random rois CV-10fold
    if strcmp(parms.skeleton.lambda,'cv')
        disp('running 20-rois CV-10 based lambda selection...')
        rng(0,'twister'); % for reproducibility during tests
        n_cases = 20;
        cases_idx = randperm(n_rois,n_cases);
        lambdas_test = zeros(1,n_cases);
        % compute lasso cv for each case
        for c = 1:n_cases
            y_id = cases_idx(c);
            y = D(:,y_id);
            y = zscore(y);
            X_ids = 1:n_rois;
            X_ids(y_id) = [];
            X = D(:,X_ids);
            X = zscore(X,0,1);
            % to run the cv in parallel. Not working in my machine.
            %opts = statset('UseParallel',true);
            [~, FitInfo] = lasso(X,y,'Intercept',false,"UseCovariance",true,'CV',10); %matlab'Options',opts);
            lambdas_test(1,c) = FitInfo.LambdaMinMSE;
        end
        % the chosen lambda is the average across the cases
        lambda = mean(lambdas_test);
    else
        % use the provided lambda
        lambda = str2double(parms.skeleton.lambda);
    end
        
    % save lambda somewhere? necessary
    disp(['Using lambda = ',num2str(lambda)]);
    
    % create a pool of workers
    % parallel workers are not working in my machine.
    % need to check, why...check in amarel
    % parpool(2)
    
    for r = 1:n_rois
        % target roi
        y_id = r;
        y = D(:,y_id);
        y = zscore(y);
        % regressors rois
        X_ids = 1:n_rois;
        X_ids(y_id) = [];
        X = D(:,X_ids);
        X = zscore(X,0,1);
        % lasso with a given lambda
        B = lasso(X,y,'Intercept',false,'UseCovariance',true,'Lambda',lambda);
        % Use a temporary variable to store the result and assign outside of parfor
        temp_fc = zeros(1, n_rois);  % Create a temporary row vector
        % Assign the computed lasso coefficients to the temporary variable
        temp_fc(X_ids) = B;  % Place the lasso results in the corresponding positions
        % Assign the result to the fc matrix
        fc(r, :) = temp_fc;  
        
    end
    warning('on','stats:lasso:StandardizeIgnored');
    
    
    % Transform FC data into a list of tuples:
    % source(src),target(trg),fc-value(fc)
    % symmetrize the fc weights by taking the average of (i,j) and (j,i)
    fc = (fc + fc')/2;
    % Extract unique values; below the diagonal
    nrChannels = size(fc,1);
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