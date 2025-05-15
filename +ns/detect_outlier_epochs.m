function flags = detect_outlier_epochs(data, srate, options)
%DETECT_OUTLIER_EPOCHS Flags EEG epochs based on multiple artifact criteria using vectorized operations and nested functions.
%
%   [bad_epochs_idx, bad_epochs_flags] = DETECT_OUTLIER_EPOCHS(data, srate, Name, Value, ...)
%   analyzes the input EEG data (Epochs x Channels x Timepoints) and flags
%   epochs containing artifacts based on specified criteria. It uses a parallel,
%   vectorized approach where possible, with nested functions for clarity.
%
%   Required Input Arguments:
%       data (Epochs x Channels x TimePoints): The preprocessed EEG data.
%       srate (double): Sampling rate of the data in Hz.
%
%   Optional Name-Value Pair Arguments:
%       'flat_threshold_sd' (double): Minimum std dev across time for any
%           channel to be considered non-flat. Default: 1 (uV).
%       'amplitude_threshold_peak' (double): Maximum peak-to-peak
%           amplitude allowed for any channel within an epoch. Default: 150 (uV).
%       'variance_z_threshold' (double): Robust Z-score threshold for
%           maximum channel variance within an epoch (flags if Z > threshold). Default: 5.
%       'hf_cutoff_hz' (double): Cutoff frequency in Hz for defining
%           high-frequency noise. Default: 50 (Hz).
%       'hf_z_threshold' (double): Robust Z-score threshold for maximum
%           channel high-frequency noise power (flags if Z > threshold). Default: 5.
%       'correlation_z_threshold' (double): Robust Z-score threshold for
%           mean absolute pairwise channel correlation (flags if Z > threshold). Default: 3.
%       'flag_if_channels_noisy' (vector): Vector of channel indices. If provided,
%           an epoch is only flagged by flat, amplitude, variance, or HF noise
%           criteria if one of THESE channels triggers the flag. Correlation
%           flagging is unaffected. Z-scoring for variance/HF is still based
%           on the max across ALL channels. Default: [].
%
%   Output Arguments:
%       bad_epochs_idx (logical vector): A logical vector (n_epochs x 1)
%           where 'true' indicates the epoch is flagged as bad.
%       bad_epochs_flags (logical matrix): A logical matrix (n_epochs x 5)
%           detailing which criterion/criteria caused each epoch rejection.
%           Columns correspond to:
%           1: Flat Signal
%           2: Amplitude Threshold
%           3: High Variance (Z-score > Threshold)
%           4: HF Noise (Z-score > Threshold)
%           5: High Correlation (Z-score > Threshold)
%
%   Example Usage:
%       % Assume 'eeg_data' is 200 epochs x 64 chans x 1000 pts
%       % Assume sampling rate is 500 Hz
%       fs = 500;
%       n_epochs = 200; n_chans = 64; n_pts = 1000;
%       eeg_data = randn(n_epochs, n_chans, n_pts); % Example data
%
%       % Flag only if channels 5 or 10 are noisy (for applicable criteria)
%       noisy_channels = [5, 10];
%       [bad_idx, details] = detect_outlier_epochs(eeg_data, fs, ...
%                               'amplitude_threshold_peak', 200, ...
%                               'flag_if_channels_noisy', noisy_channels);
%
%       fprintf('Detected %d bad epochs out of %d.\n', sum(bad_idx), size(eeg_data, 1));
%       % good_data = eeg_data(~bad_idx, :, :); % Keep only good epochs

arguments (Input)
    % Data dimensions: Epochs x Channels x Timepoints
    data (:,:,:) double {mustBeNumeric}
    srate (1,1) double {mustBeNumeric, mustBePositive} % Required positional

    % Optional Name-Value arguments collected into 'options' structure
    options.flat_threshold_sd (1,1) double {mustBeNumeric, mustBeNonnegative} = .1 % Default 1 uV
    options.amplitude_threshold_peak (1,1) double {mustBeNumeric, mustBePositive} = 250 % Default 150 uV
    options.variance_z_threshold (1,1) double {mustBeNumeric, mustBePositive} = 3.291 % Default 5 Z
    options.hf_cutoff_hz (1,1) double {mustBeNumeric, mustBePositive} = 50 % Default 50 Hz
    options.hf_z_threshold (1,1) double {mustBeNumeric, mustBePositive} = 3.291 % Default 5 Z
    options.correlation_z_threshold (1,1) double {mustBeNumeric, mustBePositive} = 3.291 % Default 3 Z
    options.criterion_channels (1,:) double {mustBeNumeric, mustBeInteger, mustBePositive} = [] % Default empty
    
    options.epoch_no (:,1) {mustBeNumeric, mustBeInteger, mustBePositive} = [] % if trial_no is different than epoch index, provide this to return correct trial numbers
    options.exclude (:,1) {mustBeNumeric, mustBeInteger, mustBePositive} = [] % Epochs to exclude from analysis
end
arguments (Output)
    flags (1, 1) ns.NSFlags
end


% --- Get Data Dimensions ---
[n_epochs, n_channels, n_timepoints] = size(data);
% exclude epochs if provided
if isempty(options.epoch_no)
    options.epoch_no = 1:n_epochs;
end

original_epoch_no = setdiff(options.epoch_no, options.exclude);

included_idx = ~ismember(options.epoch_no, options.exclude);
data = data(included_idx,:,:);
n_epochs = sum(included_idx);

fprintf('Analyzing %d epochs, %d channels, %d timepoints (Fs=%.1f Hz).\n', n_epochs, n_channels, n_timepoints, srate);

% --- Validate flag_if_channels_noisy ---
if ~isempty(options.criterion_channels)
    if any(options.criterion_channels > n_channels) || any(options.criterion_channels < 1)
        error('Values in flag_if_channels_noisy must be valid channel indices between 1 and %d.', n_channels);
    end
    fprintf('  Flagging based on specified noisy channels: %s\n', mat2str(options.criterion_channels));
else
    options.criterion_channels = 1:n_channels;
end

% --- Sanity Check Dependent Parameters ---
    if options.hf_cutoff_hz >= srate / 2
        warning('HF cutoff frequency (%.1f Hz) is >= Nyquist frequency (%.1f Hz). Check parameters.', options.hf_cutoff_hz, srate/2);
    end

% --- Initialization ---
flags = ns.NSFlags(options);

% =========================================================================
% --- Calculate Flags using Nested Functions ---
% =========================================================================
fprintf('Calculating flags...\n');

% --- 1. Flat Signal ---
isFlagged = flag_flat_signal();
flags.flat = original_epoch_no(isFlagged);
fprintf('  Flat signal check complete with %d of %d epochs flagged.\n', sum(isFlagged), n_epochs);

% --- 2. Amplitude Threshold ---
isFlagged = flag_amplitude();
flags.amplitude = original_epoch_no(isFlagged);
fprintf('  Amplitude threshold check complete with %d of %d epochs flagged.\n', sum(isFlagged), n_epochs);

% --- 3. High Variance ---
isFlagged = flag_variance();
flags.variance = original_epoch_no(isFlagged);
fprintf('  Variance check complete with %d of %d epochs flagged.\n', sum(isFlagged), n_epochs);

% --- 4. HF Noise ---
isFlagged = flag_hf_noise();
flags.HFNoise = original_epoch_no(isFlagged);
fprintf('  HF noise check complete with %d of %d epochs flagged.\n', sum(isFlagged), n_epochs);

% --- 5. High Correlation ---
isFlagged = flag_correlation();
flags.correlation = original_epoch_no(isFlagged);
fprintf('  Correlation check complete with %d of %d epochs flagged.\n', sum(isFlagged), n_epochs);

% =========================================================================
% --- Nested Functions for Criteria Evaluation ---
% =========================================================================

    % --- 1. Flag Flat Signal Epochs ---
    function flags = flag_flat_signal()
        % Calculate std dev across time (dim 3) for each epoch/channel
        channel_sds = std(data, 0, 3); % Result: n_epochs x n_channels
        % Identify which channels are flat
        is_flat_per_channel = channel_sds < options.flat_threshold_sd; % n_epochs x n_channels

        % Flag epoch if ANY channel is below threshold
        flags = any(is_flat_per_channel(:, options.criterion_channels), 2);
        
    end

    % --- 2. Flag High Amplitude Epochs ---
    function flags = flag_amplitude()
        % Calculate peak-to-peak amplitude across time (dim 3) for each epoch/channel
        peak_to_peak = max(data, [], 3) - min(data, [], 3); % Result: n_epochs x n_channels
        % Identify which channels exceed threshold
        is_high_amp_per_channel = peak_to_peak > options.amplitude_threshold_peak; % n_epochs x n_channels

        % Flag epoch if ANY channel exceeds threshold
        flags = any(is_high_amp_per_channel(:, options.criterion_channels), 2);
            
    end

    % --- 3. Flag High Variance Epochs ---
    function flags = flag_variance()
        flags = false(n_epochs, 1); % Initialize flags
        % Calculate variance across time (dim 3) for each epoch/channel
        channel_vars = var(data, 0, 3, 'omitnan'); % Result: n_epochs x n_channels

        % Find the median variance AND the index of the channel with max variance for each epoch
        
        variance_metrics = median(channel_vars(:, options.criterion_channels), 2, 'omitnan'); % Results: n_epochs x 1

        if ~all(isnan(variance_metrics))
            % Calculate Z-scores based on max variance across ALL channels
            variance_zscores = gen.robust_z(variance_metrics, 1, 'std');
            % Identify epochs exceeding the Z-score threshold
            flags = variance_zscores > options.variance_z_threshold;
            flags(isnan(variance_zscores)) = false; % Treat NaN Z-scores as non-outliers

        else
            warning('Variance metrics were all NaN. Skipping variance check.');
        end
    end

    % --- 4. Flag High HF Noise Epochs ---
    function flags = flag_hf_noise()
        flags = false(n_epochs, 1); % Initialize flags

        NFFT = 2^nextpow2(n_timepoints);
        f = srate/2*linspace(0,1,NFFT/2+1); % Frequency vector
        hf_indices = f >= options.hf_cutoff_hz; % Logical indices for HF bins

        if ~any(hf_indices)
            warning('No frequency bins >= cutoff_hz. Check parameters/srate. Skipping HF noise check.');
            return; % Return all false flags
        end

        % Apply Hann window across time dimension (3)
        hann_win_reshaped = reshape(hann(n_timepoints)', [1, 1, n_timepoints]);
        data_windowed = data .* hann_win_reshaped;

        % Perform FFT along time dimension (3)
        fft_result = fft(data_windowed, NFFT, 3); % Result: n_epochs x n_channels x NFFT

        % Calculate single-sided power spectrum
        P2 = abs(fft_result / n_timepoints);
        P1 = P2(:,:,1:NFFT/2+1);
        P1(:,:,2:end-1) = 2*P1(:,:,2:end-1);

        % Calculate mean power in HF bins for each epoch/channel
        hf_power_per_channel = mean(P1(:,:,hf_indices), 3, 'omitnan'); % Mean along frequency dim
        hf_power_per_channel = reshape(hf_power_per_channel, [n_epochs, n_channels]);

        % Find the maximum HF power AND the index of the channel with max HF power for each epoch
        hf_noise_metrics = median(hf_power_per_channel(:, options.criterion_channels), 2, 'omitnan'); % Results: n_epochs x 1

        if ~all(isnan(hf_noise_metrics))
            % Calculate Z-scores based on max HF noise across ALL channels
            hf_noise_zscores = gen.robust_z(hf_noise_metrics, 1, 'std');
             % Identify epochs exceeding the Z-score threshold
            flags = hf_noise_zscores > options.hf_z_threshold;
            flags(isnan(hf_noise_zscores)) = false; % Treat NaN Z-scores as non-outliers
            
        else
             warning('HF Noise metrics were all NaN. Skipping HF noise check.');
        end
    end

    % --- 5. Flag High Correlation Epochs ---
    function flags = flag_correlation()
        % This check is global and not affected by flag_if_channels_noisy
        flags = false(n_epochs, 1); % Initialize flags

        if n_channels <= 1
            fprintf('  Skipping correlation calculation and check (requires >1 channel).\n');
            return; % Return all false flags
        end

        fprintf('  Calculating correlation metrics (using arrayfun)...\n');
        correlation_helper = @(epoch_data_slice) calculate_single_epoch_correlation(squeeze(epoch_data_slice), options.criterion_channels);
        epoch_indices = (1:n_epochs)';
        correlation_metrics = arrayfun(@(idx) correlation_helper(data(idx, :, :)), epoch_indices);

        if any(isnan(correlation_metrics))
             warning('%d epochs resulted in NaN correlation metric (potentially due to flat channels within epoch). These will not be flagged by correlation.', sum(isnan(correlation_metrics)));
        end

        if ~all(isnan(correlation_metrics))
            correlation_zscores = gen.robust_z(correlation_metrics, 1, 'std');
             % Handle potential NaNs in Z-score output
            nan_z_corr = isnan(correlation_zscores);
            if any(nan_z_corr)
                warning('gen.robust_z returned NaN for %d correlation scores. These epochs will not be flagged by correlation.', sum(nan_z_corr));
                correlation_zscores(nan_z_corr) = -inf; % Ensure NaNs don't exceed threshold
            end
            % Flag if Z-score is GREATER than the threshold (one-sided test)
            flags = correlation_zscores > options.correlation_z_threshold;
        else
             warning('Correlation metrics were all NaN. Skipping correlation Z-score check.');
        end
    end

% =========================================================================
end % End main function detect_outlier_epochs

% =========================================================================
% --- Helper Function for Arrayfun (Correlation - Unchanged) ---
% =========================================================================
function avg_abs_correlation = calculate_single_epoch_correlation(single_epoch_data, incld_channel_idx)
    % Input: single_epoch_data (n_channels x n_timepoints)
    [n_chans, ~] = size(single_epoch_data);

    if nargin < 2
        incld_channel_idx = 1:n_chans;
    end
    if n_chans < 2
        avg_abs_correlation = NaN;
        return;
    end

    % Check for channels with zero variance before calculating correlation
    if any(var(single_epoch_data, 0, 2) < eps)
        avg_abs_correlation = NaN;
        return;
    end

    % Calculate correlation matrix (channels are rows, so transpose)
    corr_matrix = corr(single_epoch_data'); % corr works on columns

    % Extract unique off-diagonal absolute correlation values
    mask = triu(true(size(corr_matrix)), 1); % Upper triangle, excluding diagonal
    subset = zeros(size(mask));
    subset(incld_channel_idx,:) = true;
    corr_values = abs(corr_matrix(mask & subset));

    if isempty(corr_values)
        avg_abs_correlation = NaN;
    else
        avg_abs_correlation = mean(corr_values, 'omitnan'); % Mean of absolute pairwise correlations
    end
end
