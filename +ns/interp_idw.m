function [signal, options] = interp_idw(signal, coord_or_dist, channel_dim, options)
%INTERP_IDW Spatio-temporal Inverse Distance Weighting interpolation.
%
%   Interpolates 'bad' channels based on signals from 'good' channels. 
%   Dynamically handles partially missing data (NaNs) by grouping data points 
%   with identical missing data patterns to optimize performance.
%
%   Inputs:
%       signal         - N-dimensional data array (>= 2D).
%       coord_or_dist  - Matrix of Cartesian coordinates or a Distance Matrix.
%       channel_dim    - Integer. Specifies which dimension of 'signal' 
%                        corresponds to the channels. Default: 1.
%
%   Optional Name-Value Pair Arguments:
%       'good'                   - Vector of 1-based indices for permitted predictors.
%       'bad'                    - Vector of 1-based indices for strict targets.
%       'idwPower'               - Scalar > 0. The power parameter (p). Default: 2.0.
%       'min_n_predictor'        - Integer >= 1. Minimum number of good channels required. Default: 2.
%       'max_predictor_distance' - Scalar > 0. Maximum distance allowed for a predictor to be used. Default: inf.
%       'verbose'                - Logical. If true, prints progress. Default: false.
%       'distanceMatrix'         - Logical. Set to true if coord_or_dist is a distance matrix.

arguments
    signal {mustBeNumeric, mustBeReal}
    coord_or_dist (:,:) {mustBeNumeric, mustBeReal}
    channel_dim (1,1) {mustBeNumeric, mustBeInteger, mustBePositive} = 1
    
    % Optional Name-Value Pairs:
    options.good (1,:) {mustBeNumeric, mustBeInteger, mustBePositive} = []
    options.bad (1,:) {mustBeNumeric, mustBeInteger, mustBePositive} = []
    options.idwPower (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = 2.0
    options.min_n_predictor (1,1) {mustBeNumeric, mustBeInteger, mustBePositive} = 2
    options.max_predictor_distance (1,1) {mustBeNumeric, mustBeReal, mustBePositive} = inf
    options.verbose (1,1) {mustBeNumericOrLogical} = false
    options.distanceMatrix (1,1) logical = false
end

verbose = logical(options.verbose);
idwPower = options.idwPower;
min_n_predictor = options.min_n_predictor;
max_predictor_dist = options.max_predictor_distance;
epsilon = 1e-9; 

% --- 1. Validate and Generalize N-Dimensional Array ---
original_dims = size(signal);
num_dims = length(original_dims);

assert(num_dims >= 2, 'signal must have at least 2 dimensions.');
assert(channel_dim <= num_dims, ...
    'channel_dim (%d) exceeds the number of dimensions in signal (%d).', channel_dim, num_dims);

nChannels = original_dims(channel_dim);

% --- Assertions ---
assert(size(coord_or_dist, 1) == nChannels, ...
    'Number of channels in signal dimension %d (%d) must match coord_or_dist rows (%d).', ...
    channel_dim, nChannels, size(coord_or_dist, 1));
assert(~options.distanceMatrix || isscalar(unique(size(coord_or_dist))), ...
    'Distance Matrix must be square.');
if ~isempty(options.good) && ~isempty(options.bad)
    assert(isempty(intersect(options.good, options.bad)), ...
        'Good and bad channel indices cannot overlap.');
end

% --- Flatten N-Dimensional signal to 2D (Channels x Everything Else) ---
% 1. Permute the array so the channel dimension is first
perm_vec = [channel_dim, setdiff(1:num_dims, channel_dim)];
signal_perm = permute(signal, perm_vec);

% 2. Reshape into a 2D matrix
nOther = numel(signal) / nChannels;
signal_2d = reshape(signal_perm, nChannels, nOther);

% --- 2. Identify Missing Data Patterns ---
if verbose; fprintf('INFO: Scanning data for missing value patterns...\n'); end
[uniquePatterns, patternAssign] = identify_patterns(signal_2d);
numPatterns = size(uniquePatterns, 2);

if verbose
    fprintf('INFO: Found %d unique spatial pattern(s) across %d total data points.\n', numPatterns, nOther);
end

% --- 3. Process Each Pattern ---
for p = 1:numPatterns
    % Get the boolean mask for this specific pattern (1 = NaN, 0 = Valid)
    currentPatternMask = uniquePatterns(:, p); 
    pointIdx = find(patternAssign == p);
    
    % Determine Source (Predictors) and Target (Bad) channels for this pattern
    [targets, predictors] = resolve_roles(currentPatternMask, options.good, options.bad, nChannels);
    
    % --- CHECKPOINT 1: Global Pattern Check ---
    if isempty(targets)
        continue;
    elseif length(predictors) < min_n_predictor
        warning('interpolate:InsufficientPredictors', ...
            'Pattern %d: Found %d global predictor(s), but minimum required is %d. Skipping interpolation for %d target(s) across %d data point(s).', ...
            p, length(predictors), min_n_predictor, length(targets), length(pointIdx));
        continue; 
    end
    
    % Compute Distances early for thresholding
    distances = compute_distances(coord_or_dist, targets, predictors, options.distanceMatrix);
    
    % --- CHECKPOINT 2: Local Target-Specific Distance Check ---
    valid_predictor_counts = sum(distances <= max_predictor_dist, 2);
    valid_targets_mask = valid_predictor_counts >= min_n_predictor;
    
    if any(~valid_targets_mask)
        failed_targets = targets(~valid_targets_mask);
        warning('interpolate:DistantPredictors', ...
            'Pattern %d: %d target(s) (e.g., Ch %d) lack the minimum %d predictor(s) within max distance. Left as NaN.', ...
            p, length(failed_targets), failed_targets(1), min_n_predictor);
        
        % Filter out the targets that failed the distance check
        targets = targets(valid_targets_mask);
        distances = distances(valid_targets_mask, :);
    end
    
    % If all targets in this pattern failed the check, skip to next pattern
    if isempty(targets)
        continue;
    end

    if verbose && numPatterns > 1
        fprintf('  -> Pattern %d/%d: Interpolating %d target(s) using %d predictor(s) across %d data point(s).\n', ...
            p, numPatterns, length(targets), length(predictors), length(pointIdx));
    end

    % Extract data 
    predictorData = signal_2d(predictors, pointIdx);
    
    % Perform IDW (Pass max_dist to ensure distant weights are strictly 0)
    interpValues = perform_idw_core(distances, predictorData, idwPower, epsilon, class(signal_2d), max_predictor_dist);
    
    % Places the calculated interpolated values back into the 2D signal matrix.
    signal_2d(targets, pointIdx) = interpValues;
end

% --- 4. Reconstruct Original N-Dimensional Array ---
% 1. Reshape back to the permuted dimensions
permuted_dims = original_dims(perm_vec);
signal_perm_restored = reshape(signal_2d, permuted_dims);

% 2. Inverse permute back to the user's original dimensional layout
signal = ipermute(signal_perm_restored, perm_vec);

if verbose
    fprintf('INFO: Spatio-temporal IDW interpolation complete.\n');
end

end % END OF MAIN FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL/HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [uniquePatterns, patternAssign] = identify_patterns(signal_2d)
    nanMask = isnan(signal_2d);
    [uniquePatts, ~, patternAssign] = unique(nanMask', 'rows');
    uniquePatterns = uniquePatts'; 
end

function [targets, predictors] = resolve_roles(patternMask, optGood, optBad, nChannels)
    allChans = 1:nChannels;
    hasNan = find(patternMask);       
    hasData = find(~patternMask);     
    
    if isempty(optGood) && isempty(optBad)
        targets = hasNan;
        predictors = hasData;
    elseif ~isempty(optBad) && isempty(optGood)
        targets = optBad;
        predictors = intersect(setdiff(allChans, optBad), hasData);
    elseif isempty(optBad) && ~isempty(optGood)
        targets = intersect(hasNan, setdiff(allChans, optGood));
        predictors = intersect(optGood, hasData);
    else
        targets = optBad;
        predictors = intersect(optGood, hasData);
    end
end

function distances = compute_distances(coords, targets, predictors, isDistanceMatrix)
    if isDistanceMatrix
        distances = coords(targets, predictors);
    else
        targetCoords = coords(targets, :);
        predictorCoords = coords(predictors, :);
        distances = pdist2(targetCoords, predictorCoords, 'euclidean');
    end
end

function interpolatedSignal = perform_idw_core(distances, predictorData, p, epsilon, signalClass, maxDist)
    [num_targets, ~] = size(distances);
    num_points = size(predictorData, 2);
    
    interpolatedSignal = zeros(num_targets, num_points, signalClass);

    % Calculate base weights 
    weights = 1.0 ./ (max(distances, epsilon) .^ p);
    weights(~isfinite(weights)) = 0; 
    
    % --- Apply Distance Constraints ---
    weights(distances > maxDist) = 0; 

    sumWeights = sum(weights, 2);
    zeroSumMask = abs(sumWeights) < epsilon;
    validWeightsMask = ~zeroSumMask;

    % Calculate interpolation for targets with valid weights
    if any(validWeightsMask)
        weights_valid = weights(validWeightsMask, :);
        sumWeights_valid = sumWeights(validWeightsMask);
        
        % Perform weighted sum using matrix multiplication (Spatial x Temporal)
        weightedDataSum_valid = weights_valid * predictorData;
        
        % Normalize by the sum of weights
        interpValues_valid = weightedDataSum_valid ./ sumWeights_valid;
        interpolatedSignal(validWeightsMask, :) = interpValues_valid;
    end

    % Direct coincidence logic (distance < epsilon)
    [minDists, closestPredIdx] = min(distances, [], 2);
    coincidentMask = minDists < epsilon;
    if any(coincidentMask)
        for i = find(coincidentMask)'
            interpolatedSignal(i, :) = predictorData(closestPredIdx(i), :);
        end
    end
end