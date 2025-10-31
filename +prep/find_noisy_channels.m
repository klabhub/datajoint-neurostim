function BB = find_noisy_channels(signal, fs, options)
% FIND_NOISY_CHANNELS Identifies noisy EEG channels using PREP criteria.
%   
%   badByFlags = find_noisy_channels(signal, Name, Value, ...)
%
%   Usage Examples:
%     % Minimal call (requires Fs, skips RANSAC)
%     badByFlags = find_noisy_channels(eegData, 250);
%
%     % Specify Fs and Locations (runs RANSAC by default)
%     badByFlags = find_noisy_channels(eegData, 250, 'channelLocations', coordMatrix);
%
%     % Specify all, disable RANSAC
%     badByFlags = find_noisy_channels(eegData, 250, 'channelLocations', coordMatrix, 'enableRansac', false);
%
%   Inputs:
%       signal        - (nChannels x nTimepoints) EEG data matrix.
%       varargin      - Optional Name-Value pairs:
%           'channelLocations'  - (nChannels x 3) Matrix of Cartesian coordinates
%                                 [X, Y, Z]. REQUIRED to run RANSAC. Default: [].
%            See parameter definitions inside the arguments blocks for
%            additional parms and their defaults.
%
%   Outputs:
%       badByFlags - prep.badBy object containing indices of noisy channels by category
%                       (.badByNan, .badByFlat, .badByDeviation, .badByCorrelation,
%                       .badByHFNoise, .badByRansac, .all) and parameters.
%
%   Requires: Signal Processing Toolbox (filtering, corr, mad/iqr).
%   SEE ALSO prep.badBy
arguments
    % Required Positional Argument
    signal (:,:) {mustBeNumeric, mustBeNonempty}   
    fs (1,1) double {mustBeNumeric, mustBeNonnegative}

    % Optional Name-Value Pair Arguments    
    options.channelLocations (:,3) {mustBeNumeric, mustBeReal, mustBeFinite} = []
    options.fs {mustBeNumeric, mustBeScalarOrEmpty, mustBePositive} = []
    options.highPassCutoff {mustBeNumeric, mustBeScalarOrEmpty, mustBeNonnegative} = 1

    % Amplitude Cut-Off
    options.highAmplitudeCutOff (1,1) {mustBeNumeric, mustBePositive} = 500
    options.highAmplitudeMaxFrac (1,1) {mustBeNumeric, mustBePositive, mustBeLessThanOrEqual(options.highAmplitudeMaxFrac, 1)} = 0.25

    % Deviation
    options.deviationThreshold (1,1) {mustBeNumeric, mustBePositive} = 5

    % Correlation
    options.correlationFrequencyCutoff (1,1) {mustBeNumeric, mustBePositive} = 50
    options.correlationWindowSeconds (1,1) {mustBeNumeric, mustBePositive} = 1
    options.correlationThreshold (1,1) {mustBeNumeric, mustBeInRange(options.correlationThreshold, 0, 1)} = 0.4
    options.correlationMaxBadWindows (1,1) {mustBeNumeric, mustBeInRange(options.correlationMaxBadWindows, 0, 1)} = 0.01
    options.correlationUseMovingWindow (1,1) {mustBeNumericOrLogical} = false % Flag for centered moving window (step=1) vs. non-overlapping

    % HF Noise
    options.noiseFrequencyCutoff (1,1) {mustBeNumeric, mustBePositive} = 50
    options.noiseThreshold (1,1) {mustBeNumeric, mustBePositive} = 5

    % RANSAC
    options.enableRansac (1,1) {mustBeNumericOrLogical} = true
    options.ransacWindowSeconds (1,1) {mustBeNumeric, mustBePositive} = 4
    options.ransacChannelFraction (1,1) {mustBeNumeric, mustBeInRange(options.ransacChannelFraction, 0, 1, 'exclusive')} = 0.25
    options.ransacCorrelationThresholdMethod = 'absolute' % if 'robust_z'
    options.ransacCorrelationThreshold (1,1) {mustBeNumeric} = 0.75 
    options.ransacCorrelationMethod = 'Spearman' % Rank correlation ny default
    options.ransacMaxBadWindows (1,1) {mustBeNumeric, mustBeInRange(options.ransacMaxBadWindows, 0, 1)} = 0.4
    options.ransacMinimumSampleSize (1,1) {mustBeNumeric, mustBeInteger, mustBePositive} = 5 % minimum number of channels needed as predictors
    options.rngSeed {mustBeNumeric, mustBePositive} = []

    % Interpolation
    options.interpolationMethod (1,1) string {mustBeMember(options.interpolationMethod,["inverse_distance" "spherical_spline"])} = "inverse_distance"
    options.splineParams (1,3) {mustBeNumeric, mustBeVector} = [0 4 7] % For spherical_spline
    options.distanceParams = struct(idwPower = 'optimize')   % For inverse_distance interpolation

end

[nChannels, nTimepoints] = size(signal);  
signal = double(signal);
assert(isempty(options.channelLocations) || size(options.channelLocations, 1) ==nChannels,'Number of rows in ''channelLocations'' (%d) must match the number of channels in ''signal'' (%d).', size(options.channelLocations, 1), nChannels);      
assert(~options.enableRansac || ~isempty(options.channelLocations),'RANSAC needs channelLocations.\n');
switch options.interpolationMethod
    case "spherical_spline"
        options.interp_func = @prep.interpolate_by_spherical_spline;
        options.interpParams = options.splineParams;
    case "inverse_distance"
        options.interp_func = @prep.interpolate_by_inverse_distance;
        options.interpParams = namedargs2cell(options.distanceParams);
end
options.fs = fs; % To pass to sub functions.

% --- Optional High-pass Filter ---
if ~isempty(options.highPassCutoff)  && options.highPassCutoff > 0    
    fprintf('INFO: Applying temporary %.2f Hz high-pass filter...\n', options.highPassCutoff);
    try
        nyquist = fs / 2;
        if options.highPassCutoff < nyquist
            hp_cutoff_norm = options.highPassCutoff / nyquist;
            [b, a] = butter(4, hp_cutoff_norm, 'high');
            signal = filtfilt(b, a, signal')'; % Apply filter
        else
            warning('High-pass cutoff >= Nyquist. Skipping filter.');
        end
    catch ME
        warning(ME.identifier,'Could not apply high-pass filter: %s', ME.message);
    end
end

% --- Call Artifact Detection Functions Sequentially ---
funs = {@badByNanFlat, @badByAmplitude @badByDeviation @badByHfNoise @badByCorrelation,@badByRansac};
BB = prep.badBy(options); % Initalize the flags object
channels = 1:nChannels;
goodChannels = 1:nChannels; %Initially all are assumed good.
robustMean = zeros(1,nTimepoints);
for f =1:numel(funs)
    funName = func2str(funs{f});
    fprintf('Checking %s:\n ',funName)
    thisOptions = options;
    thisOptions.channelLocations = thisOptions.channelLocations(goodChannels,:);
    thisBad = funs{f}(signal(goodChannels,:,:)-robustMean,thisOptions);
    BB.(funName) = goodChannels(thisBad);
    fprintf('Found %d bad channels.\n',sum(thisBad))
    goodChannels = setdiff(channels,BB.all);
    if isempty(goodChannels)
        fprintf('No good channels left. Exiting... \n');
        break;
    end
    if options.interpolationMethod ~="" && any(thisBad)
        %Robust referencing - recalculate the global mean 
        intpSignalForMean = options.interp_func(signal, options.channelLocations,goodChannels, BB.all, options.interpParams{:});
        robustMean = median(intpSignalForMean, 1);       
    end    
end

disp(BB)

end 

function [isBad] = badByNanFlat(signal,~)
% Detects NaN and flat channels.
nanSum = sum(isnan(signal), 2);
isNan = nanSum > 0;
stdDev = std(signal, 0, 2, 'omitnan');
isZeroVariance = stdDev < eps;
isFlat = all(abs(diff(signal, 1, 2))<eps,2);
isBad = isNan| isZeroVariance | isFlat;
end % detect_bads_by_nan_flat

% --- Amplitude Cut Off ---
function [isBad] = badByAmplitude(signal,options)
signal = signal - median(signal, 1);
isHighSample = abs(signal) > options.highAmplitudeCutOff;
isBadChannel = mean(isHighSample, 2) > options.highAmplitudeMaxFrac;
isBad =isBadChannel;
end

% --- Deviation ---
function [isBad] = badByDeviation(signal,options)
% Detects channels with abnormally high amplitude/deviation.
try
    channelIQR = iqr(signal, 2);
    channelRobustStd = 0.7413 * channelIQR;
    channelRobustStd(channelRobustStd < eps) = eps;
    deviationZ = do.robust_z(channelRobustStd, 1, "std"); % Z-score across channels in subset
    isBad = deviationZ > options.deviationThreshold;
catch ME
    warning(ME.identifier, 'Deviation check failed: %s', ME.message);
    isBad = false(size(channel,1),1,1);
end
end % detect_bads_by_deviation

% --- Correlation ---
function [isBad] = badByCorrelation(signal,options)
% Detects channels poorly correlated with other channels.
[nChannels,nTimepoints] = size(signal);
try
    nyquist = options.fs / 2; lp_cutoff = options.correlationFrequencyCutoff;
    if lp_cutoff < nyquist
        [b, a] = butter(4, lp_cutoff / nyquist, 'low');
        signal = filtfilt(b, a, double(signal'))';
        % else  No filter if cutoff >= Nyquist
    end

    windowLengthSamples = round(options.correlationWindowSeconds * options.fs);
    if windowLengthSamples < 2 || windowLengthSamples > nTimepoints
        error('Correlation window size invalid for signal length.');
    end

    if options.correlationUseMovingWindow
        nWindows = nTimepoints;
        win_mask_func = @(idx) [idx - floor(windowLengthSamples/2), idx - ceil(windowLengthSamples/2)];

    else % Non-overlapping
        nWindows = floor(nTimepoints / windowLengthSamples);
        win_mask_func = @(idx) [(idx - 1) * windowLengthSamples + 1, idx*windowLengthSamples];
    end
    if nWindows < 1, error('Signal too short for correlation windows.'); end

    maxCorrelationsWindowed = zeros(nChannels, nWindows);
    for iWin = 1:nWindows

        win_idx = win_mask_func(iWin);
        winStart = win_idx(1);
        winEnd = win_idx(2);
        if winStart <= 0 || winEnd > nWindows*windowLengthSamples
            % edge handling with moving windows
            continue;
        end
        windowData = signal(:, winStart:winEnd);
        windowData(~isfinite(windowData)) = 0;
        if size(windowData,1) <= 1, corrMatrix = 1; else, corrMatrix = corrcoef(windowData'); end
        if any(isnan(corrMatrix(:))) || size(corrMatrix,1) ~= size(windowData,1), continue; end % Skip bad window

        for iChan = 1:nChannels % Index within subset
            otherChansMask = true(1, nChannels);
            otherChansMask(iChan) = false;
            if ~any(otherChansMask), maxCorrelationsWindowed(iChan, iWin) = 0; continue; end % Only 1 channel case
            corrsWithOthers = abs(corrMatrix(iChan, otherChansMask));
            if isempty(corrsWithOthers) || all(isnan(corrsWithOthers))
                maxCorrelationsWindowed(iChan, iWin) = 0;
            else
                maxCorrelationsWindowed(iChan, iWin) = prctile(corrsWithOthers(~isnan(corrsWithOthers)), 98);
            end
            if isempty(maxCorrelationsWindowed(iChan, iWin))
                maxCorrelationsWindowed(iChan, iWin) = 0;
            end % Handle prctile edge case
        end % iChan loop
    end % iWin loop



    lowCorrFraction = sum(maxCorrelationsWindowed < options.correlationThreshold, 2, 'omitnan') / nWindows;
    isBad = (lowCorrFraction > options.correlationMaxBadWindows);
catch ME
    warning(ME.identifier, 'Correlation check failed: %s', ME.message);
    isBad = false(nChannels,1,1) ;
end
end % detect_bads_by_correlation

% --- HF Noise ---
function [isBad] = badByHfNoise(signal,options)
% Detects channels with excessive high-frequency noise relative to low-freq.
try
    nyquist = options.fs / 2;
    if options.noiseFrequencyCutoff < nyquist
        [b_noise, a_noise] = butter(4, options.noiseFrequencyCutoff  / nyquist, 'low');
        lowFreqComponent = filtfilt(b_noise, a_noise, signal')';

        highFreqComponent = signal - lowFreqComponent;
        madHighFreq = mad(highFreqComponent, 1, 2); % Use function handle
        madLowFreq = mad(lowFreqComponent, 1, 2);   % Use function handle
        madLowFreq(madLowFreq < eps) = eps;
        noisinessRatio = madHighFreq ./ madLowFreq;
        if any(~isfinite(noisinessRatio))
            noisinessRatio(~isfinite(noisinessRatio)) = max(noisinessRatio(isfinite(noisinessRatio)), 1);
        end
        noiseZ = do.robust_z(noisinessRatio, 1, "std"); % Z-score across subset
        isBad = noiseZ > options.noiseThreshold;
    else
        isBad = false(size(signal,1),1);
    end
catch ME
    warning(ME.identifier, 'HF Noise check failed: %s', ME.message);
    isBad = false(size(signal,1),1);
end
end % detect_bads_by_hf_noise


% --- RANSAC ---
function isBad = badByRansac(signal,options)
% Detects channels poorly predicted by their neighbors using RANSAC.
isBad = [];
% Check enough channels remain
[nChannels,nTimepoints] =size(signal);
if nChannels <  options.ransacMinimumSampleSize + 1
    fprintf('Not enough good channels (%d) remaining for RANSAC (minimum %d required). Skipping.', nChannels, options.ransacMinimumSampleSize + 1);    
    return;
end
% Check window size validity
ransacWindowLengthSamples = round(options.ransacWindowSeconds * options.fs);
if ransacWindowLengthSamples > nTimepoints || ransacWindowLengthSamples < 2
    fprintf('RANSAC window size (%.2f s = %d samples) invalid for signal length %d. Skipping RANSAC.', ...
        options.ransacWindowSeconds, ransacWindowLengthSamples, nTimepoints);
    return;
end
nRansacWindows = floor(nTimepoints / ransacWindowLengthSamples);
if nRansacWindows < 1
    fprintf('Signal too short for RANSAC windows. Skipping RANSAC.'); 
    return; 
end


fprintf(' INFO: Using %d RANSAC windows of length %d samples (%.2f s).\n', ...
    nRansacWindows, ransacWindowLengthSamples, options.ransacWindowSeconds);

% --- RANSAC Core Logic ---
ransacBadWindowCount = zeros(nChannels, 1); % Store counts using original indices
ransacTested = false(nChannels, 1); % Keep track of channels tested
ransacCorrCoeff = nan(nChannels, nRansacWindows);
% Use original signal for RANSAC prediction targets/data TODO - Preprocessed filter? 
ransacSignal = signal;

try
    % Outer Loop: Iterate through each potentially good channel as target
    for iTargetChan = 1:nChannels        
        ransacTested(iTargetChan) = true;
        if mod(iTargetChan, 10) == 0 || iTargetChan == nChannels
            fprintf('  RANSAC: first %d target channel out of %d are complete...\n', ...
                iTargetChan-1, nChannels);
        end

        predictorPool = setdiff(1:nChannels, iTargetChan);
        % Calculate based on fraction of CURRENTLY good channels
        nPredictorsNeeded = ceil(options.ransacChannelFraction * nChannels);
        nPredictorsNeeded = min(nPredictorsNeeded, length(predictorPool)); % Ensure we don't request more than available

        if nPredictorsNeeded < options.ransacMinimumSampleSize
            fprintf('   Skipping channel %d: Not enough predictors (%d available, min %d needed based on fraction/sample size).\n', ...
                iTargetChan, length(predictorPool), options.ransacMinimumSampleSize);
            continue; % Skip to next target channel
        end

        % Inner Loop: Iterate through time windows
        for iWin = 1:nRansacWindows
            winStart = (iWin - 1) * ransacWindowLengthSamples + 1;
            winEnd = iWin * ransacWindowLengthSamples;
            windowTimeIndices = winStart:winEnd;

            % Randomly select predictor subset for this window
            predictorIdx = predictorPool(randperm(length(predictorPool), nPredictorsNeeded));

            % Get data for this window (use original signal)
            windowData = ransacSignal(:, windowTimeIndices);

            % Predict target channel using spherical spline
            try
                % Call using function handle defined in parent scope
                interpolatedSignal = options.interp_func(windowData, options.channelLocations, predictorIdx, iTargetChan, options.interpParams{:});
                interpolatedSignal = interpolatedSignal(iTargetChan,:);
            catch ME_spline
                warning('Spline interpolation failed for chan %d, win %d: %s. Counting as bad window.', iTargetChan, iWin, ME_spline.message);
                ransacBadWindowCount(iTargetChan) = ransacBadWindowCount(iTargetChan) + 1;
                continue; % Skip correlation for this window
            end

            % Calculate correlation
            actualSignal = windowData(iTargetChan, :);
            if std(actualSignal, 0, 2) < eps || std(interpolatedSignal, 0, 2) < eps % Handle flat lines
                correlation = 0;
            else
                try % Use column vectors for corr
                    correlation = corr(actualSignal(:), interpolatedSignal(:), 'Type', options.ransacCorrelationMethod);
                catch ME_corr
                    warning('Correlation calculation failed for chan %d, win %d: %s. Counting as bad window.', iTargetChan, iWin, ME_corr.message);
                    correlation = NaN;
                end
            end

            % Check threshold
            if strcmp(options.ransacCorrelationThresholdMethod, 'absolute') && (isnan(correlation) || correlation < options.ransacCorrelationThreshold)
                ransacBadWindowCount(iTargetChan) = ransacBadWindowCount(iTargetChan) + 1;
            end
            ransacCorrCoeff(iTargetChan, iWin) = correlation;
        end % End window loop
    end % End target channel loop

    % If method is robust_z, first check the distribution of
    % correlations
    if strcmp(options.ransacCorrelationThresholdMethod, 'robust_z')
        [z_corrCoeff] = do.robust_z(ransacCorrCoeff,2);
        ransacBadWindowCount = sum(z_corrCoeff < options.ransacCorrelationThreshold, 2, "omitnan");
    end

    % Classify based on fraction of bad windows for tested channels    
    isBad =  (ransacBadWindowCount/nRansacWindows > options.ransacMaxBadWindows) & ransacTested;
catch ME_ransac
    warning(ME_ransac.identifier, 'RANSAC check failed during execution: %s', ME_ransac.message);
    isBad = false(nChannels,1);
end

end % detect_bads_by_ransac
