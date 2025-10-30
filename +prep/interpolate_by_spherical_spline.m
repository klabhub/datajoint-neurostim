function interpolatedData = interpolate_by_spherical_spline(signalData, channelCoords, goodChanIdx, badChanIdx, interpParams, verbose)
% CUSTOM_SPHERICAL_SPLINE Interpolates channel data using EEGLAB's core logic.
%
%   interpolatedData = custom_spherical_spline(signalData, channelCoords, goodChanIdx, badChanIdx, interpParams)
%
%   This function extracts the core spherical spline logic from EEGLAB's
%   eeg_interp.m to interpolate specified channels based on others, without
%   requiring the EEGLAB EEG structure format.
%
%   Inputs:
%       signalData     - Data matrix (channels x time).
%       channelCoords  - Matrix of 3D Cartesian coordinates (X, Y, Z) for
%                        ALL channels (num_channels x 3). Row index must
%                        match signalData row index.
%       goodChanIdx    - Vector of indices for the 'good' channels to use
%                        for interpolation (the predictors).
%       badChanIdx     - Vector of indices for the 'bad' channels to be
%                        interpolated (the targets).
%       interpParams   - (Optional) [1x3 array] Parameters [lambda m n] for
%                        the spline calculation. Defaults to [0 4 7] matching
%                        EEGLAB's 'spherical' default (Perrin et al., 1989).
%                          lambda: Regularization parameter (often 0 or small).
%                          m:      Spline order exponent (e.g., 4).
%                          n:      Max degree for Legendre polynomials (e.g., 7).
%
%   Outputs:
%       interpolatedData - Matrix (length(badChanIdx) x time) containing the
%                          interpolated data for the bad channels.
%
%   Based on code from EEGLAB's eeg_interp.m by Arnaud Delorme. Requires the
%   local functions spheric_spline_core and computeg_core (copied below).

if nargin < 4
    error('Requires at least signalData, channelCoords, goodChanIdx, and badChanIdx.');
end
if nargin < 5 || isempty(interpParams)
    interpParams = [0 4 7]; % Default [lambda m n] for 'spherical'
    fprintf('INFO: Using default spherical spline parameters [lambda=%g, m=%d, n=%d]\n', interpParams(1), interpParams(2), interpParams(3));
end

if nargin < 6

    verbose = false;

end

if length(interpParams) ~= 3
    error('interpParams array must have 3 elements [lambda m n]');
end
if isempty(goodChanIdx)
    error('Need at least one good channel index for interpolation.');
end
if isempty(badChanIdx)
    warning('No bad channel indices provided. Returning empty matrix.');
    interpolatedData = zeros(0, size(signalData, 2), class(signalData));
    return;
end
if any(isnan(channelCoords(:)))
     error('Channel coordinates contain NaNs. Cannot proceed.');
end
if max(max(goodChanIdx), max(badChanIdx)) > size(channelCoords, 1) || min(min(goodChanIdx), min(badChanIdx)) < 1
    error('Channel indices are out of bounds for the provided coordinate matrix.');
end

rad = sqrt(sum(channelCoords.^2, 2));
if any(abs(rad) < eps), error('Good channel coordinates include origin [0,0,0].'); end
channelCoords = channelCoords./rad;
% --- Prepare Coordinates (Extract and Normalize) ---
% Coordinates for good channels (predictors)
xelec = channelCoords(goodChanIdx, 1);
yelec = channelCoords(goodChanIdx, 2);
zelec = channelCoords(goodChanIdx, 3);

% Coordinates for bad channels (targets)
xbad = channelCoords(badChanIdx, 1);
ybad = channelCoords(badChanIdx, 2);
zbad = channelCoords(badChanIdx, 3);

% --- Prepare Data ---
% Data from good channels
values = signalData(goodChanIdx,:);

% --- Perform Spherical Spline Interpolation ---
% Calls the core logic extracted from eeg_interp.m
if verbose
    fprintf('INFO: Performing spherical spline interpolation for %d channel(s)...\n', length(badChanIdx));
end
try
    % Note: The original eeg_interp spheric_spline returns xbad,ybad,zbad as first
    % 3 outputs, which we don't need here. The 4th output is the data.
    [~, ~, ~, interpolatedData] = spheric_spline_core( xelec, yelec, zelec, xbad, ybad, zbad, values, interpParams);
catch ME
     error('custom_spherical_spline:InterpolationError', ...
          'Error during core spline calculation: %s', ME.message);
end

if verbose
    fprintf('INFO: Interpolation complete.\n');
end

signalData(badChanIdx,:) = interpolatedData;
interpolatedData = signalData;

end % END OF MAIN FUNCTION custom_spherical_spline

% =========================================================================
% == COPIED AND RENAMED LOCAL FUNCTIONS FROM EEGLAB'S eeg_interp.m ==
% == (Minimal changes for standalone use, kept original logic) ==
% =========================================================================

function [xbad, ybad, zbad, allres] = spheric_spline_core( xelec, yelec, zelec, xbad, ybad, zbad, values, params)
% Core spherical spline calculation logic from eeg_interp.
% Renamed to spheric_spline_core to avoid conflicts.

    newchans = length(xbad); % Number of channels to interpolate
    numpoints = size(values,2); % Number of time points

    % Calculate G matrices using geometric relationships and Legendre polynomials
    Gelec = computeg_core(xelec,yelec,zelec,xelec,yelec,zelec,params);
    Gsph  = computeg_core(xbad,ybad,zbad,xelec,yelec,zelec,params);

    % Compute solution for parameters C (spline coefficients)
    % -------------------------------------------------------
    meanvalues = mean(values, 1); % Mean across channels for each time point
    values = values - repmat(meanvalues, [size(values,1) 1]); % Make mean zero

    values = [values; zeros(1,numpoints)]; % Add constraint sum(c_i)=0

    lambda = params(1); % Regularization parameter
    
    % Solve the linear system G*C = V for C
    % Using pseudo-inverse (pinv) for potentially better stability
    C = pinv([Gelec + eye(size(Gelec)) * lambda; ones(1, length(Gelec))]) * values;
    % Original line used pinv - might be more robust than backslash (\)
    % C = [Gelec + eye(size(Gelec)) * lambda; ones(1, length(Gelec))] \ values; % Alternative

    clear values; % Free memory
    allres = zeros(newchans, numpoints, class(C)); % Initialize output array

    % Apply results (Calculate interpolated values)
    % -------------------------------------------
    % The interpolated value at a 'bad' location j is sum(C_i * Gsph_ji) + mean
    for j = 1:size(Gsph,1) % For each bad channel location
        % Sum product of coefficients C (channels x time) and Gsph row (1 x good_channels)
        % Need to multiply C' (time x channels) by Gsph row vector transposed (channels x 1)
        % Result should be time x 1, then transpose back? No, think matrix.
        % C is (M+1) x T, Gsph is Bad x M. We need Bad x T.
        % allres(j,:) = sum(C(1:end-1,:) .* repmat(Gsph(j,:)', [1 size(C,2)])); % Original line, seems slightly off if C includes C0?
        % Let C_coeffs = C(1:end-1, :); which is M x T
        % Let Gsph_row_j = Gsph(j, :); which is 1 x M
        % We want sum over M: Gsph_row_j * C_coeffs => (1 x M) * (M x T) = (1 x T)
        allres(j,:) = Gsph(j,:) * C; % Use matrix multiplication
    end
    allres = allres + repmat(meanvalues, [size(allres,1) 1]); % Add mean back

end % END OF spheric_spline_core

% ------------------
% compute G function (copied and renamed)
% ------------------
function g = computeg_core(x,y,z,xelec,yelec,zelec, params)
% Core computeg logic from eeg_interp

    % Calculate the cosine of the angle between points based on normalized coordinates
    % cos(angle) = dot(u,v) = u_x*v_x + u_y*v_y + u_z*v_z
    % This assumes inputs x,y,z and xelec,yelec,zelec are already normalized!
    % cos_angle = zeros(length(x), length(xelec));
    % for i = 1:length(x)
    %     for j = 1:length(xelec)
    %         cos_angle(i,j) = x(i)*xelec(j) + y(i)*yelec(j) + z(i)*zelec(j);
    %     end
    % end
    

    % Assumes x,y,z and xelec,yelec,zelec are column vectors (or handles matrix inputs correctly)
    % Combine coordinates into matrices: targets (Nbad x 3) and sources (Ngood x 3)
    targetCoords = [x(:), y(:), z(:)];     % Ensure column vectors, forms Nbad x 3 matrix
    sourceCoords = [xelec(:), yelec(:), zelec(:)]; % Ensure column vectors, forms Ngood x 3 matrix
    
    % Calculate dot products using matrix multiplication
    % Resulting cos_angle will be Nbad x Ngood
    cos_angle = targetCoords * sourceCoords'; % Transpose sourceCoords

    % Clip values just in case of numerical errors slightly outside [-1, 1]
    cos_angle(cos_angle >  1.0) =  1.0;
    cos_angle(cos_angle < -1.0) = -1.0;
    % The EI calculation in original code was:
    % EI = 1 - sqrt((repmat(x(:),1,length(xelec)) - repmat(xelec,length(x(:)),1)).^2 +...
    %              (repmat(y(:),1,length(xelec)) - repmat(yelec,length(x(:)),1)).^2 +...
    %              (repmat(z(:),1,length(xelec)) - repmat(zelec,length(x(:)),1)).^2);
    % As analyzed before, cos(angle) = 1 - distance^2 / 2.
    % Let's use the direct cos(angle) calculation above which is standard.
    % The Legendre functions should take cos(angle) as input.

    g = zeros(length(x), length(xelec)); % Size: targets x sources
    m = params(2);    % Spline order exponent
    maxn = params(3); % Max Legendre polynomial degree

    for n = 1:maxn
        % Calculate P_n(cos_angle)
        % legendre(n, X) returns a matrix P, where P(1,:) is P_n(X)
        % Need to handle matrix input X carefully depending on MATLAB version/behavior.
        % Let's compute legendre for each pair type if needed, or assume vectorized works.
        try
            Pn_cos_angle = legendre(n, cos_angle);
            if size(Pn_cos_angle, 1) > 1 % If legendre returns multiple orders...
                Pn_cos_angle = Pn_cos_angle(1, :); % ...take the first row (P_n)
                 % Need to reshape back if cos_angle was matrix
                 if size(cos_angle,1) > 1 && size(cos_angle,2) > 1
                    Pn_cos_angle = reshape(Pn_cos_angle, size(cos_angle));
                 end
            end
           
            % Check dimensions (Pn should match cos_angle)
            if ~isequal(size(Pn_cos_angle), size(cos_angle))
                 error('Legendre output size mismatch.');
                 % Fallback: loop computation? (Less efficient)
                 % Pn_cos_angle = zeros(size(cos_angle));
                 % for r=1:size(cos_angle,1)
                 %    for c=1:size(cos_angle,2)
                 %       tmpL = legendre(n, cos_angle(r,c));
                 %       Pn_cos_angle(r,c) = tmpL(1);
                 %    end
                 % end
            end
            
        catch ME_leg
             error('custom_spherical_spline:LegendreError', ...
                  'Error evaluating Legendre polynomial P_%d. Input range issue? Error: %s', n, ME_leg.message);
        end

        % Apply the summation formula for g(x)
        g_term = ((2*n + 1) / (n^m * (n + 1)^m));
        g = g + g_term * Pn_cos_angle;
    end
    g = g / (4 * pi);

end % END OF computeg_core

function distMatrix = calculate_spherical_distances(chanLocs, varargin)
% CALCULATE_SPHERICAL_DISTANCES Computes pairwise spherical distances between channel locations.
%
%   distMatrix = calculate_spherical_distances(chanLocs)
%   distMatrix = calculate_spherical_distances(chanLocs, 'Radius', r)
%   distMatrix = calculate_spherical_distances(chanLocs, ..., 'UnusedParam', value)
%
%   Calculates the great-circle (spherical) distance between all pairs of
%   channel locations provided in Cartesian coordinates.
%
%   Inputs:
%       chanLocs      - (n_ch x 3) Matrix of Cartesian coordinates [X, Y, Z]
%                       for n_ch channels. Rows correspond to channels.
%       varargin      - Optional Name-Value pairs:
%           'Radius'  - (Optional) Positive scalar specifying the radius 'r'
%                       of the sphere on which the coordinates lie. If not
%                       provided, a radius of 1.0 is assumed (unit sphere),
%                       and the output distances are effectively in radians.
%                       If provided, the output distances will be in the
%                       same units as the radius (e.g., mm). Default: 1.0.
%          
%   Outputs:
%       distMatrix    - (n_ch x n_ch) Symmetric matrix where distMatrix(i, j)
%                       is the spherical distance (arc length) between
%                       channel i and channel j. The diagonal elements are 0.
%
%   Formula Used:
%       1. Normalize coordinates to unit vectors: u_i = P_i / ||P_i||
%       2. Calculate cosine of angle: cos(theta_ij) = dot(u_i, u_j)
%       3. Calculate angle: theta_ij = acos(cos(theta_ij))
%       4. Calculate spherical distance: dist_ij = r * theta_ij
%
%   Example:
%       coords = [1 0 0; 0 1 0; -1 0 0; 0 0 1]; % Example coordinates
%       % Calculate distances on unit sphere (output in radians)
%       distMatUnit = calculate_spherical_distances(coords);
%       % Calculate distances assuming sphere radius 85 mm
%       distMatMM = calculate_spherical_distances(coords, 'Radius', 85);
%
%   See also: acos, pdist2

% --- Input Parser ---
p = inputParser;
p.FunctionName = 'calculate_spherical_distances';

addRequired(p, 'chanLocs', @(x) isnumeric(x) && ismatrix(x) && size(x, 2) == 3 && ~isempty(x));
addParameter(p, 'Radius', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, chanLocs, varargin{:});

chanLocs = p.Results.chanLocs;
chanLocs = chanLocs./sum(chanLocs.^2,2);
r = 1;%p.Results.Radius;
n_ch = size(chanLocs, 1);

if r ~= 1.0
    fprintf('INFO: Calculating spherical distances assuming radius r = %g\n', r);
else
    fprintf('INFO: Calculating spherical distances assuming unit radius (r = 1.0). Distances will be in radians.\n');
end

% --- Handle Edge Cases ---
if n_ch < 2
    distMatrix = zeros(n_ch);
    if n_ch == 1, fprintf('INFO: Only one channel provided. Returning distance matrix [0].\n'); end
    if n_ch == 0, fprintf('INFO: No channels provided. Returning empty matrix.\n'); end
    return;
end

% --- Normalize Coordinates ---
% Calculate norms (magnitude) of each coordinate vector
norms = sqrt(sum(chanLocs.^2, 2));

% Check for zero vectors (coordinates at origin)
if any(norms < eps)
    zeroIdx = find(norms < eps);
    error('Channel coordinates for indices [%s] are at or very near the origin [0,0,0]. Cannot normalize.', num2str(zeroIdx'));
end

% Normalize coordinates to lie on the unit sphere
normLocs = chanLocs ./ norms; % Uses implicit expansion

% --- Calculate Cosine of Angles ---
% Dot product of normalized vectors gives cosine of angle between them
% Result is n_ch x n_ch matrix
cos_angle_matrix = normLocs * normLocs';

% --- Clip for Numerical Stability ---
% Ensure values are strictly within [-1, 1] for acos
cos_angle_matrix(cos_angle_matrix >  1.0) =  1.0;
cos_angle_matrix(cos_angle_matrix < -1.0) = -1.0;

% --- Calculate Angles ---
% Angle theta = acos(cos(theta))
theta_matrix = acos(cos_angle_matrix); % Output is in radians

% --- Calculate Spherical Distances ---
distMatrix = r * theta_matrix;

% --- Ensure Diagonal is Zero ---
% Due to potential floating point inaccuracies, explicitly set diagonal to 0
distMatrix(logical(eye(n_ch))) = 0;

fprintf('INFO: Distance matrix calculation complete.\n');

end % END OF FUNCTION calculate_spherical_distances

function scaled_signal = scale_signal_to_neighbors(signal, chanLocs, chan_to_scale_idx, n_neighbors, varargin)
% SCALE_SIGNAL_TO_NEIGHBORS Rescales specified channels to match neighbor variance.
%
%   scaled_signal = scale_signal_to_neighbors(signal, chanLocs, chan_to_scale_idx, n_neighbors)
%   scaled_signal = scale_signal_to_neighbors(..., 'ExcludeChannels', exclude_idx)
%
%   This function rescales the time series signal of specified channels
%   (target channels) so that their variance matches the average variance
%   of their nearest spatial neighbors. Channels listed in 'ExcludeChannels'
%   are not considered as neighbors.
%
%   Inputs:
%       signal             - (n_ch x n_timepoints) Data matrix. Rows are channels,
%                            columns are time points.
%       chanLocs           - (n_ch x 3) Matrix of Cartesian coordinates [X, Y, Z]
%                            for each channel. Used to find neighbors.
%       chan_to_scale_idx  - Vector of indices of the channels in 'signal'
%                            that need to be rescaled.
%       n_neighbors        - Positive integer. The number of nearest spatial
%                            neighbors (not in ExcludeChannels) to use for
%                            calculating the reference variance.
%       varargin           - Optional Name-Value pairs:
%           'ExcludeChannels' - (Optional) Vector of channel indices to
%                               exclude from being considered as neighbors.
%                               Default: [].
%
%   Outputs:
%       scaled_signal      - (n_ch x n_timepoints) Data matrix identical to 'signal'
%                            except for the rows specified in 'chan_to_scale_idx',
%                            which have been rescaled.
%
%   Requires:
%       - calculate_spherical_distances.m function on the MATLAB path.
%
%   Algorithm:
%       1. Calculate pairwise spherical distances between all channels.
%       2. For each target channel index in 'chan_to_scale_idx':
%          a. Sort channels by distance from the target channel.
%          b. Select the 'n_neighbors' closest channels that are NOT the
%             target itself and are NOT in 'ExcludeChannels'.
%          c. Calculate the variance of the target channel's signal.
%          d. Calculate the average variance of the selected neighbor channels' signals.
%          e. Compute a scaling factor = sqrt(avg_neighbor_variance / target_variance).
%          f. Multiply the target channel's signal by the scaling factor.
%
%   Example:
%       % Assume signal (128x1000), chanLocs (128x3)
%       bad_channel_indices = [10, 45];
%       known_noisy_channels = [5, 15, 46]; % Don't use these as neighbors
%       num_neighbors_for_scaling = 5;
%       scaled_signal_data = scale_signal_to_neighbors(signal, chanLocs, ...
%                                           bad_channel_indices, num_neighbors_for_scaling, ...
%                                           'ExcludeChannels', known_noisy_channels);
%
%   See also: calculate_spherical_distances, var, mean, sort, ismember, sqrt

% --- Input Validation ---
if nargin < 4
    error('Requires at least 4 input arguments: signal, chanLocs, chan_to_scale_idx, n_neighbors');
end

[n_ch, n_timepoints] = size(signal);
if size(chanLocs, 1) ~= n_ch || size(chanLocs, 2) ~= 3
    error('chanLocs must be an n_ch x 3 matrix, where n_ch matches the first dimension of signal.');
end
if any(chan_to_scale_idx < 1) || any(chan_to_scale_idx > n_ch) || ~isvector(chan_to_scale_idx)
    error('chan_to_scale_idx must be a vector containing valid channel indices (1 to %d).', n_ch);
end
if ~isscalar(n_neighbors) || n_neighbors <= 0 || mod(n_neighbors, 1) ~= 0
    error('n_neighbors must be a positive integer.');
end
if n_neighbors >= n_ch
    error('n_neighbors (%d) must be less than the total number of channels (%d).', n_neighbors, n_ch);
end
if ~exist('calculate_spherical_distances', 'file')
    error('The required helper function "calculate_spherical_distances.m" was not found on the MATLAB path.');
end

% --- Input Parser for Optional Args ---
p = inputParser;
addParameter(p, 'ExcludeChannels', [], @(x) isempty(x) || (isnumeric(x) && isvector(x) && all(x >= 1) && all(x <= n_ch) && all(mod(x,1)==0)));
parse(p, varargin{:});
exclude_idx = unique(p.Results.ExcludeChannels(:)'); % Ensure unique row vector

% --- Initialization ---
scaled_signal = signal; % Start with a copy of the original signal
tiny_variance_threshold = 1e-12; % Threshold to consider variance as effectively zero

% --- Calculate Distances (Assume Unit Sphere for neighbor ranking) ---
fprintf('Calculating distances to find neighbors...\n');
distMatrix = calculate_spherical_distances(chanLocs);
fprintf('Distance calculation complete.\n');

% --- Loop through channels to scale ---
fprintf('Scaling %d channel(s)...\n', length(chan_to_scale_idx));
for i = 1:length(chan_to_scale_idx)
    target_idx = chan_to_scale_idx(i);
    fprintf(' Processing channel index %d...\n', target_idx);

    % --- Find Neighbors (excluding target and specified exclusions) ---
    % Get distances from the target channel to all channels
    distances_from_target = distMatrix(target_idx, :);

    % Sort distances and get original indices
    [sorted_distances, sorted_indices] = sort(distances_from_target);

    % Iterate through sorted indices to find valid neighbors
    neighbor_indices = [];
    count_found = 0;
    for j = 1:length(sorted_indices) % Iterate through ALL channels sorted by distance
        current_idx = sorted_indices(j);

        % Skip if it's the target channel itself
        if current_idx == target_idx
            continue;
        end

        % Skip if it's in the exclusion list
        if ismember(current_idx, exclude_idx)
            continue;
        end

        % If valid neighbor, add to list
        neighbor_indices = [neighbor_indices, current_idx]; %#ok<AGROW>
        count_found = count_found + 1;

        % Stop when we have enough neighbors
        if count_found == n_neighbors
            break;
        end
    end

    % Check if enough valid neighbors were found
    if count_found < n_neighbors
        warning('Found only %d valid neighbors (requested %d) for channel %d after exclusions. Scaling will be based on these %d neighbors.', count_found, n_neighbors, target_idx, count_found);
        if count_found == 0
             warning('Found ZERO valid neighbors for channel %d. Skipping scaling.', target_idx);
             continue; % Skip to next channel if no neighbors found
        end
    end

    % --- Calculate Variances ---
    try
        % Variance of the target channel's signal (along time dimension 2)
        var_target = var(signal(target_idx, :), 0, 2);

        % Variance of the selected neighbor channels' signals
        vars_neighbors = var(signal(neighbor_indices, :), 0, 2);

        % Average variance of neighbors
        mean_var_neighbors = mean(vars_neighbors);

    catch ME_var
        warning('Could not calculate variance for channel %d or its neighbors. Skipping scaling. Error: %s', target_idx, ME_var.message);
        continue; % Skip to the next channel
    end

    % --- Calculate Scaling Factor ---
    scalingFactor = 1.0; % Default to 1 (no scaling)

    if var_target < tiny_variance_threshold
        warning('Target channel %d has near-zero variance (%g). Cannot calculate scaling factor. Skipping scaling.', target_idx, var_target);
        continue; % Skip scaling for this channel
    else
         ratio = max(0, mean_var_neighbors) / var_target;
         scalingFactor = sqrt(ratio);
    end

     if ~isfinite(scalingFactor) || ~isreal(scalingFactor)
         warning('Calculated scaling factor for channel %d is not finite or real (%g). Skipping scaling.', target_idx, scalingFactor);
         scalingFactor = 1.0;
         continue;
     end

    % --- Apply Scaling ---
    scaled_signal(target_idx, :) = signal(target_idx, :) * scalingFactor;
    fprintf('  Channel %d scaled by factor %.4f (Target Var: %.2g, Avg Neighbor Var: %.2g, based on %d neighbors)\n', ...
            target_idx, scalingFactor, var_target, mean_var_neighbors, count_found);

end % End loop through channels to scale

fprintf('Signal scaling complete.\n');

end % END OF FUNCTION scale_signal_to_neighbors