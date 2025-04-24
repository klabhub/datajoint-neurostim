function [signalData, options] = interpolate_by_inverse_distance(signalData, channelCoords, goodChanIdx, badChanIdx, options)
%INTERPOLATE_BY_INVERSE_DISTANCE Interpolates channel data using Inverse Distance Weighting.
%
%   signalData = interpolate_by_inverse_distance(signalData, channelCoords, ...
%                           goodChanIdx, badChanIdx, Name, Value, ...)
%
%   Interpolates the signal for 'bad' channels based on the signals from
%   'good' channels using the Inverse Distance Weighting (IDW) algorithm,
%   and updates these channels in the returned signalData matrix. Optional
%   parameters are specified using Name-Value pairs.
%
%   Inputs:
%       signalData     - Data matrix (channels x time).
%       channelCoords  - Matrix of 2D or 3D Cartesian coordinates (X, Y, [Z])
%                        for ALL channels (num_channels x num_dims). Row index
%                        must match signalData row index.
%                        - If distanceMetric is 'euclidean', 2D or 3D is supported.
%                        - If distanceMetric is 'geodesic', 3D is required.
%       goodChanIdx    - Vector of 1-based indices for the 'good' channels to use
%                        for interpolation (the predictors).
%       badChanIdx     - Vector of 1-based indices for the 'bad' channels to be
%                        interpolated (the targets).
%
%   Optional Name-Value Pair Arguments:
%       'distanceMetric' - String specifying the distance metric:
%                          "euclidean" (default): Straight-line distance in 2D/3D space. Uses pdist2.
%                          "geodesic": Shortest path distance on the surface of the
%                                      3D convex hull of ALL channel locations. Requires 3D
%                                      coordinates. Uses internal 'calculate_surface_distances'.
%       'idwPower'     - Scalar > 0 OR the string "optimize". The power
%                        parameter (p) for IDW. If "optimize", the optimal
%                        power is found using cross-validation via fminbnd (using
%                        Euclidean distances for the optimization process itself).
%                        Default: 2.0.
%       'idwOptimizationParameters' - Struct containing parameters for the
%                        optimization process, used only if idwPower is "optimize".
%                        Fields and (Defaults):
%                           .frac_predictor_ch (0.25) - Fraction of good channels
%                                                       for CV subset.
%                           .p_min             (1.0)  - Min exponent for search.
%                           .p_max             (5.0)  - Max exponent for search.
%                        Default: struct with default field values. Note: The internal
%                                 'verbose' setting of this struct is ignored; use the
%                                 main 'verbose' argument.
%       'verbose'      - Logical. If true, prints informational messages
%                        about the interpolation process (including start/end
%                        of optimization if applicable and geodesic calculation steps).
%                        Default: false.
%
%   Outputs:
%       signalData     - The input signalData matrix with rows corresponding
%                        to badChanIdx updated with interpolated values. If
%                        badChanIdx is empty, the original signalData is returned.
%
%   Algorithm Notes:
%       - Uses epsilon = 1e-9 internally for numerical stability checks.
%       - If idwPower is "optimize", calls optimize_idwPower_via_fminbnd. The
%         optimization *always* uses Euclidean distances internally for evaluating
%         candidate 'p' values, regardless of the main 'distanceMetric' setting.
%       - If distanceMetric is "geodesic", the convex hull and surface graph
%         are computed using *all* points in channelCoords to ensure graph
%         accuracy, even if some points are not in goodChanIdx or badChanIdx.
%         The required distances between bad and good channels are then extracted.
%         This requires 3D coordinates and may fail if points are collinear/coplanar.
%
%   Example:
%       % ... setup signalData, channelCoords (3D), good_indices, bad_indices ...
%
%       % Interpolate using default Euclidean distance
%       signalData_euc = interpolate_by_inverse_distance(signalData, channelCoords, ...
%                                         good_indices, bad_indices, 'verbose', true);
%
%       % Interpolate using Geodesic distance on the hull (requires 3D coords)
%       signalData_geo = interpolate_by_inverse_distance(signalData, channelCoords, ...
%                                         good_indices, bad_indices, ...
%                                         'distanceMetric', "geodesic", 'verbose', true);
%
%       % Optimize idwPower (using Euclidean distance internally for optimization),
%       % then apply using Geodesic distance for the final interpolation
%       signalData_opt_geo = interpolate_by_inverse_distance(signalData, channelCoords, ...
%                                         good_indices, bad_indices, ...
%                                         'idwPower', "optimize", ...
%                                         'distanceMetric', "geodesic", 'verbose', true);
%
%   See also: pdist2, convhull, graph, distances, perform_idw_core, optimize_idwPower_via_fminbnd

arguments % Start arguments block
    signalData (:,:) {mustBeNumeric, mustBeReal}
    channelCoords (:,:) {mustBeNumeric, mustBeReal}
    goodChanIdx (1,:) {mustBeNumeric, mustBeInteger, mustBePositive, mustBeVector}
    badChanIdx (1,:) {mustBeNumeric, mustBeInteger, mustBePositive, mustBeVector}
    % Optional Name-Value Pairs:
    options.distanceMetric (1,1) string {mustBeMember(options.distanceMetric, ["euclidean", "geodesic"])} = "euclidean"
    % Allow idwPower to be numeric OR string initially, validate further below
    options.idwPower {mustBeNumericOrStringScalar(options.idwPower)} = 2.0
    options.idwOptimizationParameters (1,1) struct = struct() % Default empty struct
    options.verbose (1,1) {mustBeNumericOrLogical} = false % Default
end % End arguments block

% --- Extract options for use ---
idwPowerInput = options.idwPower; % Keep original input separate
verbose = logical(options.verbose); % Ensure logical type (used throughout)
distanceMetric = options.distanceMetric; % Get selected metric

% --- Set Defaults and Validate Optimization Parameters ---
% Pass the main verbose flag to the default setter ONLY if the user didn't supply
% an explicit 'verbose' field within the struct (which we now ignore anyway).
% This maintains the original behavior of the helper struct function if needed elsewhere.
optParams = setDefaultOptimizationParams(options.idwOptimizationParameters, verbose);

% --- Define internal constant ---
epsilon = 1e-9; % Small value for numerical stability

% --- Additional Cross-Argument Validations ---
[numChannels, numTimepoints] = size(signalData);
if size(channelCoords, 1) ~= numChannels
    error('Number of rows in signalData (%d) must match number of rows in channelCoords (%d).', numChannels, size(channelCoords, 1));
end
% Validate coordinate dimensions based on distance metric
if strcmpi(distanceMetric, "geodesic")
    if size(channelCoords, 2) ~= 3
         error('Geodesic distance requires 3D channel coordinates (X, Y, Z). Input has %d dimensions.', size(channelCoords, 2));
    end
elseif size(channelCoords, 2) < 2
    error('channelCoords must have at least 2 dimensions (X, Y).');
end

if any(goodChanIdx > numChannels)
    error('goodChanIdx contains indices out of bounds (1 to %d).', numChannels);
end
if ~isempty(badChanIdx) && any(badChanIdx > numChannels)
    error('badChanIdx contains indices out of bounds (1 to %d).', numChannels);
end
if isempty(goodChanIdx)
    error('Need at least one good channel index for interpolation.');
end

% Check for NaNs only in relevant coordinates
if any(isnan(channelCoords(goodChanIdx,:)), 'all')
    error('Coordinates for goodChanIdx contain NaNs. Cannot proceed.');
end
if ~isempty(badChanIdx) && any(isnan(channelCoords(badChanIdx,:)), 'all')
    error('Coordinates for badChanIdx contain NaNs. Cannot proceed.');
end

if ~isempty(intersect(goodChanIdx, badChanIdx))
    warning('interpolate_by_inverse_distance:Overlap', 'Good and bad channel index lists overlap. Ensure this is intended.');
end

% --- Handle Empty Bad Channels Case ---
if isempty(badChanIdx)
    if verbose
        fprintf('INFO: No bad channel indices provided. Returning original signalData.\n');
    end
    return; % Return unmodified signalData
end

% --- Determine IDW Power (Optimize if requested) ---
if (ischar(idwPowerInput) || isstring(idwPowerInput))
    if strcmpi(idwPowerInput, "optimize")
        if length(goodChanIdx) < 2
            error('Cannot optimize idwPower: Need at least 2 good channels for cross-validation.');
        end
        if verbose
            fprintf('INFO: Optimizing IDW power (using Euclidean distances for evaluation)...\n');
            fprintf('      Optimization Params: frac=%.2f, p_min=%.2f, p_max=%.2f\n', ...
                optParams.frac_predictor_ch, optParams.p_min, optParams.p_max);
        end
        % Call optimizer - PASS MAIN VERBOSE FLAG
        % Optimization evaluation *always* uses Euclidean distance internally
        [optimalExp, ~] = optimize_idwPower_via_fminbnd(signalData, channelCoords, goodChanIdx, ...
            optParams.frac_predictor_ch, optParams.p_min, optParams.p_max, ...
            'verbose', verbose); % Pass main verbose flag
        idwPower = optimalExp; % Use the optimized value
        options.idwPoer = idwPower;
        if verbose
            fprintf('INFO: Optimization complete. Using optimal IDW power p = %.4f for final interpolation.\n', idwPower);
        end
    else
        error('Invalid string value for idwPower: "%s". Must be "optimize" or a positive scalar.', string(idwPowerInput));
    end
elseif isnumeric(idwPowerInput) && isscalar(idwPowerInput) && idwPowerInput > 0
    idwPower = idwPowerInput; % Use the user-provided numeric value
else
    error('Invalid value for idwPower. Must be a positive scalar or the string "optimize".');
end

% --- Now idwPower holds the final numeric value to use ---

% --- Proceed with Interpolation ---
numBadChans = length(badChanIdx);
numGoodChans = length(goodChanIdx);
signalClass = class(signalData);

if verbose
    fprintf('INFO: Performing Inverse Distance Weighting (power=%.4f, eps=%.2g) for %d channel(s)...\n', idwPower, epsilon, numBadChans);
    fprintf('      Using %d good channels as predictors.\n', numGoodChans);
    fprintf('      Distance Metric: %s\n', distanceMetric);
end

% Prepare Coordinates (Needed for pdist2 if Euclidean)
goodCoords = channelCoords(goodChanIdx, :);
badCoords = channelCoords(badChanIdx, :);
goodData = signalData(goodChanIdx, :); % Data from predictors

% --- Calculate Pairwise Distances Based on Metric ---
if strcmpi(distanceMetric, "geodesic")
    % Geodesic distance calculation requires 3D check already done above
    if verbose
         fprintf('INFO: Calculating geodesic (surface) distances using convex hull...\n');
         fprintf('      Using all %d points to build the graph...\n', numChannels);
    end
    % Call local function using ALL coordinates, passing the main verbose flag
    [fullSurfaceDistances, ~, ~, ~, ~] = calculate_surface_distances(channelCoords, verbose);

    % Extract the required sub-matrix of distances between bad and good channels
    distances = fullSurfaceDistances(badChanIdx, goodChanIdx);

    % Check for Inf values indicating disconnected points in the required subset
    inf_mask = isinf(distances);
    if any(inf_mask, 'all')
         num_inf = sum(inf_mask(:));
         warning('interpolate_by_inverse_distance:GeodesicInf', ...
                 'Geodesic distance calculation resulted in %d Inf value(s) between bad/good channel pairs. Check hull integrity or point locations. These pairs will have zero weight.', num_inf);
         % perform_idw_core handles Inf distance by assigning zero weight.
    end
     if verbose
         fprintf('INFO: Geodesic distance calculation complete.\n');
     end
else % Default or explicit "euclidean"
     if verbose
         fprintf('INFO: Calculating Euclidean distances...\n');
     end
     distances = pdist2(badCoords, goodCoords, 'euclidean');
      if verbose
         fprintf('INFO: Euclidean distance calculation complete.\n');
     end
end

% Initialize Temporary Output Array
interpolatedData = zeros(numBadChans, numTimepoints, signalClass);

% Identify Coincident Channels (Distance ~ 0 based on calculated metric)
[minDists, closestGoodIndRel] = min(distances, [], 2); % Find minimum distance for each bad channel
coincidentBadMask = minDists < epsilon;
interpNeededBadMask = ~coincidentBadMask;

% Handle Coincident Channels Directly
if any(coincidentBadMask)
    coincidentBadIndicesRel = find(coincidentBadMask);
    % Find the absolute index of the closest good channel for each coincident bad channel
    closestGoodIndicesAbs = goodChanIdx(closestGoodIndRel(coincidentBadMask));
    num_coincident = length(coincidentBadIndicesRel);
    if verbose
        fprintf('  INFO: %d bad channel(s) coincide with good channel(s) (distance < %.1e). Using good channel data directly.\n', num_coincident, epsilon);
    end
    for i = 1:num_coincident
        rel_idx = coincidentBadIndicesRel(i); % Relative index within the 'bad' set
        closest_good_abs = closestGoodIndicesAbs(i); % Absolute index in signalData
        interpolatedData(rel_idx, :) = signalData(closest_good_abs, :);
        % Optional: More detailed verbose message per channel
        % if verbose
        %     fprintf('      Bad channel index %d (relative #%d) coincides with good channel index %d.\n', badChanIdx(rel_idx), rel_idx, closest_good_abs);
        % end
    end
end

% Handle Channels Requiring Interpolation
if any(interpNeededBadMask)
    num_to_interp = sum(interpNeededBadMask);
    if verbose
        fprintf('  INFO: Interpolating %d non-coincident channel(s) using weighted average...\n', num_to_interp);
    end
    interpNeededBadIndicesRel = find(interpNeededBadMask);
    interpNeededBadIndicesAbs = badChanIdx(interpNeededBadIndicesRel); % Absolute indices for potential warnings
    distancesToInterp = distances(interpNeededBadMask, :); % Distances for only the channels needing interpolation

    % Call the core calculation function, passing final idwPower and defined epsilon
    interpValues = perform_idw_core(distancesToInterp, goodData, idwPower, epsilon, signalClass, interpNeededBadIndicesAbs);

    % Assign results to the correct rows in the temporary output matrix
    interpolatedData(interpNeededBadMask, :) = interpValues;
end

% Update the original signalData matrix
signalData(badChanIdx, :) = interpolatedData;

if verbose
    fprintf('INFO: Inverse Distance Weighting interpolation complete.\n');
end

end % END OF MAIN FUNCTION interpolate_by_inverse_distance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED/LOCAL FUNCTIONS START HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geodesic Distance Calculation (Local Function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [surfaceDistances, surfaceDistances_HullOnly, hullVerticesIdx, G, K] = calculate_surface_distances(points, verbose)
%CALCULATE_SURFACE_DISTANCES Computes shortest distances between points on a 3D convex hull surface. (Local Function)
%   (Full help text omitted for brevity - see previous versions if needed)
%   Receives 'verbose' flag directly from the main function.

% --- Input Validation and Defaults ---
% Basic checks (already done in main function, but good practice)
if nargin < 2 % verbose is required from main func
     error('calculate_surface_distances: Internal error - verbose flag not passed.');
end
if ~isnumeric(points) || ~ismatrix(points) || size(points, 2) ~= 3
    error('calculate_surface_distances: Internal error - Input ''points'' must be an Nx3 numeric matrix.');
end
if size(points, 1) < 4
     warning('calculate_surface_distances:SmallInput', 'Input has fewer than 4 points. Convex hull requires at least 4 non-coplanar points.');
end
% Ensure verbose is logical
verbose = logical(verbose);

% Helper function for logging (respects passed verbose flag)
log = @(msg) fprintf(msg); % Simple logger
if ~verbose
    log = @(msg) []; % No-op if not verbose
end

numVertices = size(points, 1);
% Log message adjusted slightly assuming it's called internally
log(sprintf('      [Geodesic] Input: %d points.\n', numVertices));

% --- 1. Compute Convex Hull and Triangulation ---
log('      [Geodesic] Computing convex hull...');
try
    qhull_options = {};
    if ~verbose % Suppress qhull output if main function is not verbose
        qhull_options = {'Quiet'};
    end
    [K, ~] = convhull(points(:,1), points(:,2), points(:,3), qhull_options{:});
    log(sprintf(' Done (%d faces).\n', size(K, 1)));
catch ME
    if strcmp(ME.identifier, 'MATLAB:qhullmx:DegenerateData') || contains(ME.message, 'collinear') || contains(ME.message, 'coplanar')
        error('calculate_surface_distances:DegenerateData', 'Input points appear to be collinear or coplanar. Cannot form a valid 3D convex hull surface for distance calculation.');
    else
        rethrow(ME);
    end
end

if size(K,1) < 2 % Check for meaningful hull
    warning('calculate_surface_distances:DegenerateHull', 'Convex hull is degenerate (likely flat or linear). Geodesic distances may not be meaningful.');
    surfaceDistances = inf(numVertices, numVertices);
    surfaceDistances(1:numVertices+1:end) = 0;
    hullVerticesIdx = unique(K(:));
    surfaceDistances_HullOnly = surfaceDistances(hullVerticesIdx, hullVerticesIdx);
    G = graph();
    return;
end

% --- 2. Build the Graph Representation ---
log('      [Geodesic] Building graph representation...');
adjMatrix = inf(numVertices, numVertices);
adjMatrix(1:numVertices+1:end) = 0;
edge_count = 0;
for i = 1:size(K, 1)
    face_indices = K(i, :);
    pairs = [face_indices(1), face_indices(2); face_indices(2), face_indices(3); face_indices(3), face_indices(1)];
    for p = 1:size(pairs, 1)
        idx1 = pairs(p, 1); idx2 = pairs(p, 2);
        if isinf(adjMatrix(idx1, idx2))
            point1 = points(idx1, :); point2 = points(idx2, :);
            edge_length = norm(point1 - point2);
            adjMatrix(idx1, idx2) = edge_length;
            adjMatrix(idx2, idx1) = edge_length;
            edge_count = edge_count + 1;
        end
    end
end
log(sprintf(' Done (%d unique surface edges).\n', edge_count));

% --- 3. Create Graph Object ---
log('      [Geodesic] Creating graph object...');
[sparse_i, sparse_j, sparse_v] = find(adjMatrix);
valid_edges = isfinite(sparse_v);
if ~any(valid_edges)
     warning('calculate_surface_distances:NoEdges', 'Could not find any valid edges for the graph.');
     G = graph();
else
    G = graph(sparse_i(valid_edges), sparse_j(valid_edges), sparse_v(valid_edges), numVertices);
end
log(' Done.\n');

% --- 4. Compute All-Pairs Shortest Paths ---
if isempty(G.Edges)
     warning('calculate_surface_distances:EmptyGraph', 'Resulting graph object has no edges.');
     surfaceDistances = inf(numVertices, numVertices);
     surfaceDistances(1:numVertices+1:end) = 0;
     hullVerticesIdx = unique(K(:));
     surfaceDistances_HullOnly = surfaceDistances(hullVerticesIdx, hullVerticesIdx);
else
    log('      [Geodesic] Calculating all-pairs shortest paths...');
    surfaceDistances = distances(G);
    log(' Done.\n');
    log('      [Geodesic] Extracting hull-specific results...');
    hullVerticesIdx = unique(K(:));
    surfaceDistances_HullOnly = surfaceDistances(hullVerticesIdx, hullVerticesIdx);
    log(' Done.\n');
    % Optional: Check for disconnected components
    if any(isinf(surfaceDistances_HullOnly(:))) && ~isempty(surfaceDistances_HullOnly) % Added ~isempty check
        warning('calculate_surface_distances:DisconnectedGraph', 'Some hull vertices are disconnected in the surface graph (Inf distances found). Check points/hull.');
    end
end
log('      [Geodesic] Function calculate_surface_distances finished.\n');

end % END OF LOCAL FUNCTION calculate_surface_distances

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function to set default optimization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optParamsOut = setDefaultOptimizationParams(optParamsIn, mainVerbose)
% Sets defaults for optimization parameters.
% Now IGNORES 'verbose' field within optParamsIn, uses mainVerbose if needed.
arguments
    optParamsIn (1,1) struct
    mainVerbose (1,1) logical % Required argument now
end

defaults = struct(...
    'frac_predictor_ch', 0.25, ...
    'p_min',             1.0,  ...
    'p_max',             5.0  ... % Removed 'verbose' from defaults
    );
optParamsOut = defaults;

% Check if input struct is not just the default empty one
if numel(fieldnames(optParamsIn)) > 0 % Check if user provided any fields
    userFields = fieldnames(optParamsIn);
    for i = 1:length(userFields)
        fieldName = userFields{i};
        if isfield(defaults, fieldName) % Check against known, non-verbose fields
            if isnumeric(optParamsIn.(fieldName)) && isscalar(optParamsIn.(fieldName))
                optParamsOut.(fieldName) = optParamsIn.(fieldName);
            else
                warning('setDefaultOptimizationParams: Ignoring non-scalar numeric value for field "%s". Using default.', fieldName);
            end
        elseif strcmpi(fieldName, 'verbose')
             % Explicitly ignore verbose field within the struct
             % warning('setDefaultOptimizationParams: Ignoring "verbose" field within idwOptimizationParameters struct. Use the main function''s verbose flag.');
        else
            warning('setDefaultOptimizationParams: Ignoring unknown field "%s" in idwOptimizationParameters.', fieldName);
        end
    end
end

% Final validation on parameter values (p_min, p_max, frac)
if optParamsOut.p_max <= optParamsOut.p_min
    warning('setDefaultOptimizationParams: p_max (%.2f) <= p_min (%.2f). Resetting p_max to p_min + 4.0.', optParamsOut.p_max, optParamsOut.p_min);
    optParamsOut.p_max = optParamsOut.p_min + 4.0;
end
if optParamsOut.p_min <= 0
    warning('setDefaultOptimizationParams: p_min (%.2f) must be positive. Resetting to default 1.0.', optParamsOut.p_min);
    optParamsOut.p_min = 1.0;
    if optParamsOut.p_max <= optParamsOut.p_min
        optParamsOut.p_max = optParamsOut.p_min + 4.0;
    end
end
if optParamsOut.frac_predictor_ch <= 0 || optParamsOut.frac_predictor_ch > 1
    warning('setDefaultOptimizationParams: frac_predictor_ch (%.2f) invalid. Resetting to default 0.25.', optParamsOut.frac_predictor_ch);
    optParamsOut.frac_predictor_ch = 0.25;
end
end % END OF LOCAL FUNCTION setDefaultOptimizationParams

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function for arguments block validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mustBeNumericOrStringScalar(arg)
% Validate that input is numeric scalar or char/string scalar
if ~( (isnumeric(arg) && isscalar(arg)) || ...
        (ischar(arg) && (isrow(arg) || isempty(arg))) || ...
        (isstring(arg) && isscalar(arg)) )
    error('Value must be a numeric scalar, a character row vector, or a string scalar.');
end
end % END OF LOCAL FUNCTION mustBeNumericOrStringScalar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization Function (Receives main verbose flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [optimalExp, finalCost] = optimize_idwPower_via_fminbnd(signalData, channelCoords, goodChanIdx, frac_predictor_ch, p_min, p_max, options)
%OPTIMIZE_IDWPOWER_VIA_FMINBND Finds optimal IDW power using fminbnd. (Local Function)
%   Uses Euclidean distances internally for evaluation.
arguments
    signalData (:,:) {mustBeNumeric, mustBeReal}
    channelCoords (:,:) {mustBeNumeric, mustBeReal}
    goodChanIdx (1,:) {mustBeNumeric, mustBeInteger, mustBePositive, mustBeVector}
    frac_predictor_ch (1,1) {mustBeNumeric, mustBeReal, mustBeInRange(frac_predictor_ch, 0, 1, 'exclude-lower')}
    p_min (1,1) {mustBeNumeric, mustBeReal, mustBePositive}
    p_max (1,1) {mustBeNumeric, mustBeReal} % Validated below
    % Optional Name-Value Pair:
    options.verbose (1,1) {mustBeNumericOrLogical} = false % Receives main verbose flag
end

% --- Extract options ---
verbose = logical(options.verbose); % Use the passed-in verbose flag

% --- Define internal constant ---
epsilon = 1e-9;

% --- Additional Cross-Argument Validations ---
[numChannels, ~] = size(signalData);
if size(channelCoords, 1) ~= numChannels
    error('optimize_idwPower_via_fminbnd: Signal/Coord row mismatch.');
end
if any(goodChanIdx > numChannels)
    error('optimize_idwPower_via_fminbnd: goodChanIdx out of bounds.');
end
if p_max <= p_min
    error('optimize_idwPower_via_fminbnd: p_max (%.2f) must be greater than p_min (%.2f).', p_max, p_min);
end
if length(goodChanIdx) < 2
    error('optimize_idwPower_via_fminbnd: Need at least 2 good channels.');
end
if any(isnan(channelCoords(goodChanIdx,:)), 'all')
    error('optimize_idwPower_via_fminbnd: Coordinates for goodChanIdx contain NaNs.');
end

% --- Set Optimization Options ---
if verbose
    displayMode = 'iter'; % Show fminbnd iterations if main verbose is true
else
    displayMode = 'off';
end
options_fminbnd = optimset('Display', displayMode, 'TolX', 1e-4);

% --- Create Objective Function Handle ---
% Pass epsilon explicitly. evaluate_idwPower uses Euclidean distance internally.
objectiveHandle = @(p) evaluate_idwPower(p, signalData, channelCoords, goodChanIdx, frac_predictor_ch, epsilon);

if verbose
    fprintf('INFO: Starting optimization using fminbnd in range [%.2f, %.2f]...\n', p_min, p_max);
end

% --- Run Optimization ---
[optimalExp, finalCost] = fminbnd(objectiveHandle, p_min, p_max, options_fminbnd);

% No need for separate completion message here, main function handles it.

end % END OF LOCAL FUNCTION optimize_idwPower_via_fminbnd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Core IDW Calculation (Local Function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function interpolatedSignal = perform_idw_core(distances, predictorData, p, epsilon, signalClass, targetIdsForWarning)
%PERFORM_IDW_CORE Calculates the weighted sum for IDW interpolation. (Local Function)
arguments
    distances (:,:) {mustBeNumeric, mustBeReal, mustBeNonnegative}
    predictorData (:,:) {mustBeNumeric, mustBeReal}
    p (1,1) {mustBeNumeric, mustBeReal}
    epsilon (1,1) {mustBeNumeric, mustBeReal, mustBePositive}
    signalClass (1,:) char {mustBeMember(signalClass, {'double', 'single'})}
    targetIdsForWarning (1,:) {mustBeNumeric, mustBeInteger, mustBePositive, mustBeVector} = [] % Optional absolute indices for warning
end
% --- Core Logic ---
[num_targets, num_predictors] = size(distances);
[num_predictors_data, num_time] = size(predictorData);
if num_predictors ~= num_predictors_data
    error('perform_idw_core: Mismatch between distance dimensions (%d predictors) and predictorData dimensions (%d predictors).', num_predictors, num_predictors_data);
end

interpolatedSignal = zeros(num_targets, num_time, signalClass);

% Calculate weights (handle potential Inf distances from geodesic by resulting in 0 weight)
weights = 1.0 ./ (max(distances, epsilon) .^ p);
weights(~isfinite(weights)) = 0; % Set non-finite weights (e.g., from Inf distance) to 0

sumWeights = sum(weights, 2);

% Identify targets where sum of weights is effectively zero
zeroSumMask = abs(sumWeights) < epsilon;
validWeightsMask = ~zeroSumMask;

% Calculate interpolation for targets with valid weights
if any(validWeightsMask)
    weights_valid = weights(validWeightsMask, :);
    sumWeights_valid = sumWeights(validWeightsMask);

    % Perform weighted sum using matrix multiplication
    weightedDataSum_valid = weights_valid * predictorData;

    % Normalize by the sum of weights
    interpValues_valid = weightedDataSum_valid ./ sumWeights_valid;

    interpolatedSignal(validWeightsMask, :) = interpValues_valid;
end

% Handle targets where the sum of weights was zero
if any(zeroSumMask)
    problematicTargetIndicesRel = find(zeroSumMask);
    num_problematic = length(problematicTargetIndicesRel);
    % Display absolute indices if provided
    if ~isempty(targetIdsForWarning) && length(targetIdsForWarning) == num_targets
        problematicTargetIndicesAbs = targetIdsForWarning(problematicTargetIndicesRel);
        warning('perform_idw_core: Sum of weights near-zero for target indices: [%s] (p=%.2f). All predictors were likely too distant or Inf. Interpolated values set to zero.', num2str(problematicTargetIndicesAbs(:)'), p);
    else % Fallback if absolute indices weren't passed correctly
        warning('perform_idw_core: Sum of weights near-zero for %d target channel(s) (p=%.2f). All predictors were likely too distant or Inf. Interpolated values set to zero.', num_problematic, p);
    end
    % Values already initialized to zero, so no action needed
end

end % END OF LOCAL FUNCTION perform_idw_core

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective Function for Optimization (Local Function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = evaluate_idwPower(p, signalData, channelCoords, allGoodChanIdx, frac_predictor_ch, epsilon)
%EVALUATE_IDWPOWER Objective function for IDW power optimization. (Local Function)
%   Uses Euclidean distance (pdist2) for cross-validation splits.
arguments
    p (1,1) {mustBeNumeric, mustBeReal}
    signalData (:,:) {mustBeNumeric, mustBeReal}
    channelCoords (:,:) {mustBeNumeric, mustBeReal} % Can be 2D or 3D
    allGoodChanIdx (1,:) {mustBeNumeric, mustBeInteger, mustBePositive, mustBeVector}
    frac_predictor_ch (1,1) {mustBeNumeric, mustBeReal}
    epsilon (1,1) {mustBeNumeric, mustBeReal, mustBePositive}
end
% --- Core Logic ---
numTotalGood = length(allGoodChanIdx);
n_subset_pool = max(2, round(numTotalGood * frac_predictor_ch));
n_subset_pool = min(n_subset_pool, numTotalGood); % Cannot exceed total good

if n_subset_pool < 2 % Need at least one predictor and one target
    warning('evaluate_idwPower: Too few channels in subset pool (%d) to perform split. Returning max cost.', n_subset_pool);
    cost = 2.0; return;
end

% Randomly select a subset of good channels for this fold
subset_indices = randsample(allGoodChanIdx, n_subset_pool, false);

% Randomly select one channel from the subset to be the target
target_idx_rel = randi(length(subset_indices)); % Index within subset_indices
target_channel_abs = subset_indices(target_idx_rel); % Absolute index in signalData

% The rest of the subset are predictors
predictor_set_abs = subset_indices;
predictor_set_abs(target_idx_rel) = []; % Remove the target

if isempty(predictor_set_abs)
    % This should only happen if n_subset_pool was 1, which is caught above,
    % but added as a safeguard.
    warning('evaluate_idwPower: Predictor set is empty after selecting target (subset size %d). Returning max cost.', n_subset_pool);
    cost = 2.0; return;
end

% Get coordinates and data for this fold
targetCoords = channelCoords(target_channel_abs, :);
predictorCoords = channelCoords(predictor_set_abs, :);
predictorData = signalData(predictor_set_abs, :);
actualTargetSignal = signalData(target_channel_abs, :); % Ground truth

% --- Calculate *Euclidean* distances for this CV fold ---
distances = pdist2(targetCoords, predictorCoords, 'euclidean');

minDist = min(distances); % Check if target coincides with a predictor in this fold

if minDist < epsilon
    % If coincident, use the closest predictor's data directly
    [~, closestPredIndRel] = min(distances); % Index within predictor_set_abs
    interpolatedSignal = predictorData(closestPredIndRel, :);
else
    % Otherwise, perform IDW using the current power 'p'
    interpolatedSignal = perform_idw_core(distances, predictorData, p, epsilon, class(signalData), target_channel_abs);
end

% Calculate cost (1 - correlation) if signals are not constant
stdTarget = std(actualTargetSignal);
stdInterp = std(interpolatedSignal);

if stdTarget > epsilon && stdInterp > epsilon % Use epsilon instead of 1e-12 for consistency
    corrMatrix = corrcoef(actualTargetSignal, interpolatedSignal);
    correlation = corrMatrix(1, 2);
    correlation = max(-1.0, min(1.0, correlation)); % Clamp correlation
    cost = 1.0 - correlation; % Cost is 1 minus correlation
else
    % warning('evaluate_idwPower: Target %d actual or interpolated signal is constant for p=%.2f. Cannot calculate correlation. Returning max cost.', target_channel_abs, p);
    cost = 2.0; % Assign high cost if correlation cannot be computed
end

% Handle potential NaN cost
if isnan(cost)
    % warning('evaluate_idwPower: Cost calculation resulted in NaN for target %d, p=%.2f. Returning max cost.', target_channel_abs, p);
    cost = 2.0; % Assign high cost if NaN occurs
end

end % END OF LOCAL FUNCTION evaluate_idwPower