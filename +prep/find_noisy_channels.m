function noisyChannels = find_noisy_channels(signal, options)
% FIND_NOISY_CHANNELS Identifies noisy EEG channels using PREP criteria.
%   Refactored version with nested functions and optional inputs.
%
%   noisyChannels = find_noisy_channels(signal, Name, Value, ...)
%
%   Usage Examples:
%     % Minimal call (requires Fs, skips RANSAC)
%     noisyInfo = find_noisy_channels(eegData, 'Fs', 250);
%
%     % Specify Timepoints (Fs calculated), skip RANSAC
%     noisyInfo = find_noisy_channels(eegData, 'Timepoints', timeVec);
%
%     % Specify Fs and Locations (runs RANSAC by default)
%     noisyInfo = find_noisy_channels(eegData, 'Fs', 250, 'ChannelLocations', coordMatrix);
%
%     % Specify all, disable RANSAC
%     noisyInfo = find_noisy_channels(eegData, 'Fs', 250, 'ChannelLocations', coordMatrix, 'EnableRansac', false);
%
%   Inputs:
%       signal        - (nChannels x nTimepoints) EEG data matrix.
%       varargin      - Optional Name-Value pairs:
%           'Timepoints'        - (1 x nTimepoints) Time vector. If provided and 'Fs'
%                                 is not, Fs will be calculated from this. Default: [].
%           'ChannelLocations'  - (nChannels x 3) Matrix of Cartesian coordinates
%                                 [X, Y, Z]. REQUIRED to run RANSAC. Default: [].
%           'Fs'                - Sampling frequency (Hz). REQUIRED if 'Timepoints'
%                                 is not provided or unusable for calculation. Default: [].
%           'highPassCutoff'    - Optional temporary high-pass filter cutoff (Hz). Default: 1.
%           --- Thresholds & Parameters ---
%           'deviationThreshold', 'correlation*', 'noise*', 'EnableRansac', 'ransac*'
%           (See parameter definitions inside inputParser for defaults)
%
%   Outputs:
%       noisyChannels - Structure containing indices of noisy channels by category
%                       (.badByNan, .badByFlat, .badByDeviation, .badByCorrelation,
%                       .badByHFNoise, .badByRansac, .badByLowSNR, .all) and parameters.
%
%   Requires: Signal Processing Toolbox (filtering, corr, mad/iqr).
%             ns.interpolate_by_spherical_spline.m function (if running RANSAC).
%             A robust z-score function (e.g., calculate_robust_z defined below or do.robust_z).
%             A MAD function (e.g., built-in mad or custom_mad defined below).
% --- Argument Validation using arguments block ---
arguments
    % Required Positional Argument
    signal (:,:) {mustBeNumeric, mustBeNonempty}

    % Optional Name-Value Pair Arguments
    options.Timepoints = []
    options.ChannelLocations (:,3) {mustBeNumeric, mustBeReal, mustBeFinite} = []
    options.Fs {mustBeNumeric, mustBeScalarOrEmpty, mustBePositive} = []
    options.highPassCutoff {mustBeNumeric, mustBeScalarOrEmpty, mustBeNonnegative} = 1

    % Amplitude Cut-Off
    options.highAmplitudeCutOff (1,1) {mustBeNumeric, mustBePositive} = 500
    options.highAmplitudeMaxRatio (1,1) {mustBeNumeric, mustBePositive, mustBeLessThanOrEqual(options.highAmplitudeMaxRatio, 1)} = 0.25

    % Deviation
    options.deviationThreshold (1,1) {mustBeNumeric, mustBePositive} = 3.291

    % Correlation
    options.correlationFrequencyCutoff (1,1) {mustBeNumeric, mustBePositive} = 50
    options.correlationWindowSeconds (1,1) {mustBeNumeric, mustBePositive} = 1
    options.correlationThreshold (1,1) {mustBeNumeric, mustBeInRange(options.correlationThreshold, 0, 1)} = 0.4
    options.correlationMaxBadWindows (1,1) {mustBeNumeric, mustBeInRange(options.correlationMaxBadWindows, 0, 1)} = 0.1
    options.correlationUseMovingWindow (1,1) {mustBeNumericOrLogical} = false % Flag for centered moving window (step=1) vs. non-overlapping

    % HF Noise
    options.noiseFrequencyCutoff (1,1) {mustBeNumeric, mustBePositive} = 50
    options.noiseThreshold (1,1) {mustBeNumeric, mustBePositive} = 3.291

    % RANSAC
    options.EnableRansac (1,1) {mustBeNumericOrLogical} = true
    options.ransacWindowSeconds (1,1) {mustBeNumeric, mustBePositive} = 4
    options.ransacChannelFraction (1,1) {mustBeNumeric, mustBeInRange(options.ransacChannelFraction, 0, 1, 'exclusive')} = 0.1
    options.ransacCorrelationThresholdMethod = 'absolute' % if 'robust_z'
    options.ransacCorrelationThreshold (1,1) {mustBeNumeric} = 0.5 %, mustBeInRange(options.ransacCorrelationThreshold, 0, 1)} = 0.5
    options.ransacCorrelationMethod = 'Spearman' % Rank correlation ny default
    options.ransacMaxBadWindows (1,1) {mustBeNumeric, mustBeInRange(options.ransacMaxBadWindows, 0, 1)} = 0.5
    options.ransacMinimumSampleSize (1,1) {mustBeNumeric, mustBeInteger, mustBePositive} = 5 % minimum number of channels needed as predictors
    options.rngSeed {mustBeNumeric, mustBePositive} = []

    % Interpolation
    options.interpolationMethod (1,1) string {mustBeMember(options.interpolationMethod,["inverse_distance" "spherical_spline"])} = "inverse_distance"
    options.splineParams (1,3) {mustBeNumeric, mustBeVector} = [0 4 7] % For spherical_spline
    options.distanceParams = struct(idwPower = 'optimize')   % For inverse_distance interpolation

end
% --- End of Argument Validation ---

[nChannels, nTimepoints] = size(signal); % Use dimensions from validated signal

% Perform validations that depend on the size of 'signal'
assert(isempty(options.ChannelLocations) || size(options.ChannelLocations, 1) ==nChannels,'Number of rows in ''ChannelLocations'' (%d) must match the number of channels in ''signal'' (%d).', size(options.ChannelLocations, 1), nChannels);

% --- Initial Setup ---
noisyChannels = struct();
% Handle Timepoints and Sampling Rate (Fs)
timepoints_in = options.Timepoints; % Might be empty
if isempty(options.Fs)
    if ~isempty(timepoints_in) && length(timepoints_in) >= 2
        Fs = 1 / median(diff(timepoints_in));
        if isnan(Fs) || Fs <= 0
            error('Could not calculate valid Fs from provided Timepoints. Please provide Fs directly via the ''Fs'' parameter.');
        end
        fprintf('INFO: Sampling frequency calculated as %.2f Hz from Timepoints.\n', Fs);
    else
        error('Fs not specified and cannot be calculated from Timepoints (either not provided or too short). Please provide Fs directly via the ''Fs'' parameter.');
    end
else
    Fs = options.Fs;
    fprintf('INFO: Using provided sampling frequency %.2f Hz.\n', Fs);
end
noisyChannels.Fs = Fs;

% Handle Channel Locations and RANSAC flag
channelCoords = options.ChannelLocations; % Might be empty
runRansac = options.EnableRansac && ~isempty(channelCoords); % Check if possible AND enabled
if options.EnableRansac && isempty(channelCoords)
    fprintf('INFO: RANSAC enabled but ChannelLocations not provided. Skipping RANSAC check.\n');
elseif ~options.EnableRansac
    fprintf('INFO: RANSAC check disabled by user.\n');
end

% --- Select/Check Required External/Helper Functions ---
if exist('+gen/robust_z','file')
    robust_z_func = @do.robust_z; fprintf('INFO: Using do.robust_z.\n');
elseif exist('calculate_robust_z', 'file') % Check for helper defined below
    robust_z_func = @calculate_robust_z; fprintf('INFO: Using calculate_robust_z.\n');
else
    error('Robust Z-score function not found.');
end

if exist('mad','file')
    mad_func = @mad; fprintf('INFO: Using built-in mad.\n');
else
    error('MAD function not found.');
end

interp_func = options.interpolationMethod; % Assign function handle if needed
switch options.interpolationMethod
    case "spherical_spline"
        interp_func = @ns.prep.interpolate_by_spherical_spline;
        interpParams = options.splineParams;
    case "inverse_distance"
        interp_func = @ns.prep.interpolate_by_inverse_distance;
        interpParams = namedargs2cell(options.distanceParams);
end


% --- Optional High-pass Filter ---
processedSignal = double(signal); % Use double for processing
if ~isempty(options.highPassCutoff) && options.highPassCutoff > 0
    fprintf('INFO: Applying temporary %.2f Hz high-pass filter...\n', options.highPassCutoff);
    try
        nyquist = Fs / 2;
        if options.highPassCutoff < nyquist
            hp_cutoff_norm = options.highPassCutoff / nyquist;
            [b, a] = butter(4, hp_cutoff_norm, 'high');
            processedSignal = filtfilt(b, a, processedSignal')'; % Apply filter
        else
            warning('High-pass cutoff >= Nyquist. Skipping filter.');
        end
    catch ME
        warning(ME.identifier,'Could not apply high-pass filter: %s', ME.message);
    end
end

% --- Call Detection Functions Sequentially ---

[noisyChannels.badByNan, noisyChannels.badByFlat] = detect_bads_by_nan_flat();
noisyInterim = unique([noisyChannels.badByNan, noisyChannels.badByFlat]);
goodChannelsForAmp_Idx = setdiff(1:nChannels, noisyInterim);

interim_mean = median(processedSignal(goodChannelsForAmp_Idx,:), 1);
processedSignal = processedSignal - interim_mean;
[noisyChannels.badByAmplitude, ~] = detect_bad_by_amplitude(goodChannelsForAmp_Idx);
noisyInterim = unique([noisyInterim, noisyChannels.badByAmplitude]);
goodChannelsForDeviation_Idx = setdiff(1:nChannels, noisyInterim);

processedSignal = processedSignal + interim_mean;
if ~isempty(noisyInterim)
    intpSignalForMean = interp_func(processedSignal, channelCoords, ...
        goodChannelsForDeviation_Idx, noisyInterim, interpParams{:});
    interim_mean = median(intpSignalForMean, 1);
end
processedSignal = processedSignal - interim_mean;
[noisyChannels.badByDeviation, ~] = detect_bads_by_deviation(goodChannelsForDeviation_Idx);
noisyInterim = unique([noisyInterim, noisyChannels.badByDeviation]);
goodChannelsForHFNoise_Idx = setdiff(1:nChannels, noisyInterim);

processedSignal = processedSignal + interim_mean;
if ~isempty(noisyInterim)
    intpSignalForMean = interp_func(processedSignal, channelCoords, ...
        goodChannelsForHFNoise_Idx, noisyInterim, interpParams{:});
    interim_mean = median(intpSignalForMean, 1);
end
processedSignal = processedSignal - interim_mean;
[noisyChannels.badByHFNoise, ~] = detect_bads_by_hf_noise(goodChannelsForHFNoise_Idx);
noisyInterim = unique([noisyInterim, noisyChannels.badByDeviation, ...
    noisyChannels.badByHFNoise]);
goodChannelsForCorrelation_Idx = setdiff(1:nChannels, noisyInterim);

processedSignal = processedSignal + interim_mean;
if ~isempty(noisyInterim)
    intpSignalForMean = interp_func(processedSignal, channelCoords, ...
        goodChannelsForCorrelation_Idx, noisyInterim, interpParams{:});
    interim_mean = median(intpSignalForMean, 1);
end
processedSignal = processedSignal - interim_mean;
[noisyChannels.badByCorrelation, ~] = detect_bads_by_correlation(goodChannelsForCorrelation_Idx);
noisyInterim = unique([noisyInterim, noisyChannels.badByDeviation, ...
    noisyChannels.badByCorrelation, noisyChannels.badByHFNoise]);
goodChannelsForRansac_Idx = setdiff(1:nChannels, noisyInterim);

if options.EnableRansac

    if isempty(options.rngSeed)
        options.rngSeed = randi(2^32-1);
    end
    rng(options.rngSeed);
    processedSignal = processedSignal + interim_mean;

    if ~isempty(noisyInterim)
        [intpSignalForMean, paramsN] = interp_func(processedSignal, channelCoords, ...
            goodChannelsForRansac_Idx, noisyInterim, interpParams{:});
        interim_mean = median(intpSignalForMean, 1);
    end
    prev_idwPower_idx = find(cellfun(@(x) (isstring(x) || ischar(x)) && strcmp(x,"idwPower"), interpParams)) + 1;
    if strcmp(options.interpolationMethod, "inverse_distance") && isstring(interpParams{prev_idwPower_idx}) && strcmp(interpParams{prev_idwPower_idx}, "optimize")

        % use the idwPower which was estimated at the intrep call
        interpParams{prev_idwPower_idx} = paramsN.idwPower;
        options.ransacDistanceParams.actualIDWPower = paramsN.idwPower;
    end
    processedSignal = processedSignal - interim_mean;
    noisyChannels.badByRansac = detect_bads_by_ransac(goodChannelsForRansac_Idx);

else
    noisyChannels.badByRansac = [];
end
% --- Combine Results ---
noisyChannels.badByLowSNR = intersect(noisyChannels.badByCorrelation, noisyChannels.badByHFNoise);
noisyChannels.all = unique([noisyChannels.badByNan, noisyChannels.badByFlat, ...
    noisyChannels.badByAmplitude, noisyChannels.badByDeviation, noisyChannels.badByCorrelation, ...
    noisyChannels.badByRansac, noisyChannels.badByHFNoise])';

noisyChannels.parameters = options; % Store parameters used
% --- Summary ---
fprintf('--- Summary ---\n');
fprintf(' Total unique noisy channels identified: %d\n', length(noisyChannels.all));
fprintf(' Indices: %s\n', mat2str(noisyChannels.all'));
fprintf(' Breakdown:\n');
fprintf('  NaN: %d\n', length(noisyChannels.badByNan));
fprintf('  Flat: %d\n', length(noisyChannels.badByFlat));
fprintf('  Deviation: %d\n', length(noisyChannels.badByDeviation));
fprintf('  Correlation: %d\n', length(noisyChannels.badByCorrelation));
fprintf('  HF Noise: %d\n', length(noisyChannels.badByHFNoise));
fprintf('  RANSAC: %d\n', length(noisyChannels.badByRansac));
fprintf('  Low SNR (also in Corr/HF): %d\n', length(noisyChannels.badByLowSNR));
fprintf('---------------\n');


%% --- Nested Function Definitions ---
% These functions access variables from the main function's scope:
% processedSignal, Fs, nChannels, nTimepoints, params, channelCoords,
% robust_z_func, mad_func, spline_func, runRansacCheck

    function [badNan_Idx, badFlat_Idx] = detect_bads_by_nan_flat()
        % Detects NaN and flat channels.
        fprintf('Checking for NaN/Flat channels...\n');
        nanSum = sum(isnan(processedSignal), 2);
        badNan_Idx = find(nanSum > 0)';
        stdDev = std(processedSignal, 0, 2, 'omitnan');
        isZeroVariance = stdDev < 1e-10; % Check overall variance first
        % Check for flat sections only if overall variance isn't already zero
        isFlatSection = false(nChannels, 1);
        nonZeroVarChans = find(~isZeroVariance);
        if ~isempty(nonZeroVarChans)
            diffSignal = diff(processedSignal(nonZeroVarChans,:), 1, 2);
            isFlatSection(nonZeroVarChans) = any(abs(diffSignal) < 1e-10, 2);
        end
        badFlat_all = find(isFlatSection | isZeroVariance);
        badFlat_all = badFlat_all(:)'; % make row
        badFlat_Idx = setdiff(badFlat_all, badNan_Idx); % Exclude NaNs
        fprintf(' Found %d NaN, %d Flat channels.\n', length(badNan_Idx), length(badFlat_Idx));
    end % detect_bads_by_nan_flat

% --- Amplitude Cut Off ---
    function [badAmp_Idx, mask] = detect_bad_by_amplitude(currentGood_Idx)

        fprintf('Checking amplitude criterion...\n');
        if isempty(currentGood_Idx), fprintf(' INFO: Skipping amplitude check (no good channels).\n'); return; end
        originalIndicesSubset = currentGood_Idx;
        mask = false(size(processedSignal));
        mask(currentGood_Idx,:) = abs(processedSignal(currentGood_Idx,:)) > options.highAmplitudeCutOff;
        badAmp_Idx = originalIndicesSubset((mean(mask(currentGood_Idx,:), 2) > options.highAmplitudeMaxRatio));
        badAmp_Idx = badAmp_Idx(:)'; % make row

        fprintf(' Found %d channels failing amplitude criterion.\n', length(badAmp_Idx));
    end

% --- Deviation ---
    function [badDeviation_Idx, mask] = detect_bads_by_deviation(currentGood_Idx)
        % Detects channels with abnormally high amplitude/deviation.
        fprintf('Checking deviation criterion...\n');
        badDeviation_Idx = [];
        if isempty(currentGood_Idx), fprintf(' INFO: Skipping deviation (no good channels).\n'); return; end

        signalSubset = processedSignal;
        signalSubset(currentGood_Idx,:) = NaN;
        originalIndicesSubset = currentGood_Idx;
        try
            channelIQR = iqr(signalSubset, 2);
            channelRobustStd = 0.7413 * channelIQR;
            channelRobustStd(channelRobustStd < eps) = eps;
            deviationZ = robust_z_func(channelRobustStd, 1, "std"); % Z-score across channels in subset
            mask = deviationZ > options.deviationThreshold;
            badDeviation_Idx = originalIndicesSubset(mask); % Map back
            badDeviation_Idx = badDeviation_Idx(:)'; % make row
        catch ME
            warning(ME.identifier, 'Deviation check failed: %s', ME.message);
        end
        fprintf(' Found %d channels failing deviation criterion.\n', length(badDeviation_Idx));
    end % detect_bads_by_deviation

% --- Correlation ---
    function [badCorrelation_Idx, mask] = detect_bads_by_correlation(currentGood_Idx)
        % Detects channels poorly correlated with other channels.
        fprintf('Checking correlation criterion...\n');
        badCorrelation_Idx = [];
        if isempty(currentGood_Idx), fprintf(' INFO: Skipping correlation (no good channels).\n'); return; end
        signalSubset = processedSignal(currentGood_Idx, :);
        originalIndicesSubset = currentGood_Idx;
        try
            nyquist = Fs / 2; lp_cutoff = options.correlationFrequencyCutoff;
            if lp_cutoff < nyquist
                [b, a] = butter(4, lp_cutoff / nyquist, 'low');
                filteredSignalCorr = filtfilt(b, a, double(signalSubset'))';
            else
                filteredSignalCorr = signalSubset;
            end % No filter if cutoff >= Nyquist

            windowLengthSamples = round(options.correlationWindowSeconds * Fs);
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

            maxCorrelationsWindowed = zeros(length(currentGood_Idx), nWindows);
            for iWin = 1:nWindows

                win_idx = win_mask_func(iWin);
                winStart = win_idx(1);
                winEnd = win_idx(2);
                if winStart <= 0 || winEnd > nWindows*windowLengthSamples
                    % edge handling with moving windows
                    continue;
                end
                windowData = filteredSignalCorr(:, winStart:winEnd);
                windowData(~isfinite(windowData)) = 0;
                if size(windowData,1) <= 1, corrMatrix = 1; else, corrMatrix = corrcoef(windowData'); end
                if any(isnan(corrMatrix(:))) || size(corrMatrix,1) ~= size(windowData,1), continue; end % Skip bad window

                for iChan = 1:length(currentGood_Idx) % Index within subset
                    otherChansMask = true(1, length(currentGood_Idx)); otherChansMask(iChan) = false;
                    if ~any(otherChansMask), maxCorrelationsWindowed(iChan, iWin) = 0; continue; end % Only 1 channel case
                    corrsWithOthers = abs(corrMatrix(iChan, otherChansMask));
                    if isempty(corrsWithOthers) || all(isnan(corrsWithOthers)), maxCorrelationsWindowed(iChan, iWin) = 0;
                    else, maxCorrelationsWindowed(iChan, iWin) = prctile(corrsWithOthers(~isnan(corrsWithOthers)), 98); end
                    if isempty(maxCorrelationsWindowed(iChan, iWin)), maxCorrelationsWindowed(iChan, iWin) = 0; end % Handle prctile edge case
                end % iChan loop
            end % iWin loop

            mask = nan(nChannels, nWindows);
            % mask(currentGood_Idx, :) = maxCorrelationsWindowed < options.correlationThreshold;
            lowCorrFraction = sum(maxCorrelationsWindowed < options.correlationThreshold, 2, 'omitnan') / nWindows;
            badCorrelation_Idx = originalIndicesSubset(lowCorrFraction > options.correlationMaxBadWindows); % Map back
            badCorrelation_Idx = badCorrelation_Idx(:)'; %make row

        catch ME
            warning(ME.identifier, 'Correlation check failed: %s', ME.message);
        end
        fprintf(' Found %d channels failing correlation criterion.\n', length(badCorrelation_Idx));
    end % detect_bads_by_correlation

% --- HF Noise ---
    function [badHFNoise_Idx, mask] = detect_bads_by_hf_noise(currentGood_Idx)
        % Detects channels with excessive high-frequency noise relative to low-freq.
        fprintf('Checking noisiness criterion (HF noise)...\n');
        badHFNoise_Idx = [];
        if isempty(currentGood_Idx), fprintf(' INFO: Skipping HF Noise (no good channels).\n'); return; end
        signalSubset = processedSignal(currentGood_Idx, :);
        originalIndicesSubset = currentGood_Idx;
        try
            nyquist = Fs / 2; lp_cutoff_noise = options.noiseFrequencyCutoff;
            if lp_cutoff_noise < nyquist
                [b_noise, a_noise] = butter(4, lp_cutoff_noise / nyquist, 'low');
                lowFreqComponent = filtfilt(b_noise, a_noise, double(signalSubset'))';
            else
                lowFreqComponent = signalSubset;
            end

            highFreqComponent = signalSubset - lowFreqComponent;
            madHighFreq = mad_func(highFreqComponent, 1, 2); % Use function handle
            madLowFreq = mad_func(lowFreqComponent, 1, 2);   % Use function handle
            madLowFreq(madLowFreq < eps) = eps;
            noisinessRatio = madHighFreq ./ madLowFreq;
            if any(~isfinite(noisinessRatio))
                noisinessRatio(~isfinite(noisinessRatio)) = max(noisinessRatio(isfinite(noisinessRatio)), 1);
            end
            noiseZ = robust_z_func(noisinessRatio, 1, "std"); % Z-score across subset
            mask = zeros(nChannels, nTimepoints);
            badHFNoise_Idx = originalIndicesSubset(noiseZ > options.noiseThreshold); % Map back
            badHFNoise_Idx = badHFNoise_Idx(:)'; % make row
        catch ME
            warning(ME.identifier, 'HF Noise check failed: %s', ME.message);
        end
        fprintf(' Found %d channels failing HF noise criterion.\n', length(badHFNoise_Idx));
    end % detect_bads_by_hf_noise

% --- RANSAC ---
    function badRansac_Idx = detect_bads_by_ransac(currentGood_Idx)
        % Detects channels poorly predicted by their neighbors using RANSAC.
        badRansac_Idx = [];
        % Check flag set in main function scope (based on param AND channelCoords)
        if ~runRansac
            fprintf('INFO: Skipping RANSAC check (disabled or channel locations not provided).\n');
            return;
        end
        fprintf('Checking predictability criterion (RANSAC)...\n');

        % Check enough channels remain
        minChannelsForRansac = options.ransacMinimumSampleSize + 1; % Need predictors + target
        if isempty(currentGood_Idx) || length(currentGood_Idx) < minChannelsForRansac
            warning('Not enough good channels (%d) remaining for RANSAC (minimum %d required). Skipping.', length(currentGood_Idx), minChannelsForRansac);
            return;
        end

        % Check window size validity
        ransacWindowLengthSamples = round(options.ransacWindowSeconds * Fs);
        if ransacWindowLengthSamples > nTimepoints || ransacWindowLengthSamples < 2
            warning('RANSAC window size (%.2f s = %d samples) invalid for signal length %d. Skipping RANSAC.', ...
                options.ransacWindowSeconds, ransacWindowLengthSamples, nTimepoints);
            return;
        end
        nRansacWindows = floor(nTimepoints / ransacWindowLengthSamples);
        if nRansacWindows < 1, warning('Signal too short for RANSAC windows. Skipping RANSAC.'); return; end
        fprintf(' INFO: Using %d RANSAC windows of length %d samples (%.2f s).\n', ...
            nRansacWindows, ransacWindowLengthSamples, options.ransacWindowSeconds);

        % --- RANSAC Core Logic ---
        ransacBadWindowCount = zeros(nChannels, 1); % Store counts using original indices
        ransacTested = false(nChannels, 1); % Keep track of channels tested
        ransacCorrCoeff = nan(nChannels, nRansacWindows);
        % Use original signal for RANSAC prediction targets/data
        ransacSignal = signal;

        try
            % Outer Loop: Iterate through each potentially good channel as target
            for iTargetChan_local = 1:length(currentGood_Idx)
                targetChanIdx = currentGood_Idx(iTargetChan_local); % Get original index
                ransacTested(targetChanIdx) = true;
                if mod(iTargetChan_local, 10) == 0 || iTargetChan_local == length(currentGood_Idx)
                    fprintf('  RANSAC: first %d target channel out of %d are complete. Continuing with Channel E%d)...\n', ...
                        iTargetChan_local-1, length(currentGood_Idx), targetChanIdx);
                end

                predictorPool = setdiff(currentGood_Idx, targetChanIdx);
                % Calculate based on fraction of CURRENTLY good channels
                nPredictorsNeeded = ceil(options.ransacChannelFraction * length(currentGood_Idx));
                nPredictorsNeeded = min(nPredictorsNeeded, length(predictorPool)); % Ensure we don't request more than available

                if nPredictorsNeeded < options.ransacMinimumSampleSize
                    fprintf('   Skipping channel %d: Not enough predictors (%d available, min %d needed based on fraction/sample size).\n', ...
                        targetChanIdx, length(predictorPool), options.ransacMinimumSampleSize);
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
                        interpolatedSignal = interp_func(double(windowData), channelCoords, predictorIdx, targetChanIdx, interpParams{:});
                        interpolatedSignal = interpolatedSignal(targetChanIdx,:);
                    catch ME_spline
                        warning('Spline interpolation failed for chan %d, win %d: %s. Counting as bad window.', targetChanIdx, iWin, ME_spline.message);
                        ransacBadWindowCount(targetChanIdx) = ransacBadWindowCount(targetChanIdx) + 1;
                        continue; % Skip correlation for this window
                    end

                    % Calculate correlation
                    actualSignal = windowData(targetChanIdx, :);
                    if std(actualSignal, 0, 2) < eps || std(interpolatedSignal, 0, 2) < eps % Handle flat lines
                        correlation = 0;
                    else
                        try % Use column vectors for corr
                            correlation = corr(actualSignal(:), interpolatedSignal(:), 'Type', options.ransacCorrelationMethod);
                        catch ME_corr
                            warning('Correlation calculation failed for chan %d, win %d: %s. Counting as bad window.', targetChanIdx, iWin, ME_corr.message);
                            correlation = NaN;
                        end
                    end

                    % Check threshold

                    if strcmp(options.ransacCorrelationThresholdMethod, 'absolute') && (isnan(correlation) || correlation < options.ransacCorrelationThreshold)
                        ransacBadWindowCount(targetChanIdx) = ransacBadWindowCount(targetChanIdx) + 1;
                    end
                    ransacCorrCoeff(targetChanIdx, iWin) = correlation;
                end % End window loop
            end % End target channel loop

            % If method is robust_z, first check the distribution of
            % correlations
            if strcmp(options.ransacCorrelationThresholdMethod, 'robust_z')
                % MOz: BK removed the computation of an absolute threshold
                % as it does not seem to be used (and var/med are not
                % returned by robust_z_func). Please check.
                [z_corrCoeff] = robust_z_func(ransacCorrCoeff,2);
                ransacBadWindowCount = sum(z_corrCoeff < options.ransacCorrelationThreshold, 2, "omitnan");
            end

            % Classify based on fraction of bad windows for tested channels
            %MOz : BK vectorized this. Please check.
            badRansac_Idx = find( (ransacBadWindowCount/nRansacWindows > options.ransacMaxBadWindows) & ransacTested)';
        catch ME_ransac
            warning(ME_ransac.identifier, 'RANSAC check failed during execution: %s', ME_ransac.message);
            badRansac_Idx = []; % Ensure it's empty on error
        end
        fprintf(' Found %d channels failing RANSAC criterion.\n', length(badRansac_Idx));
    end % detect_bads_by_ransac

end % END OF MAIN FUNCTION find_noisy_channels

function z = calculate_robust_z(x, dim, ~)
% Calculates robust z-score: (x - median) / (0.7413 * IQR)z

if isempty(x)
    z = [];
    return;
end

if nargin < 2 || isempty(dim)
    if isvector(x)
        dim = find(size(x) > 1);
        if isempty(dim)
            dim = 1; % Handle scalar case
        end
    else
        error('calculate_robust_z:DimensionRequired', 'Dimension must be specified for matrix input.');
    end
end

% Calculate along the specified dimension
med = median(x, dim, 'omitnan');
iqr_val = iqr(x, dim);

% Handle zero IQR to avoid division by zero
iqr_val(iqr_val < eps) = eps; % Use a small epsilon

robust_sd = 0.7413 * iqr_val;

% Use bsxfun/implicit expansion for subtraction
z = (x - med) ./ robust_sd;

% Handle cases where input was all NaN/Inf
z(~isfinite(x)) = NaN;
end

