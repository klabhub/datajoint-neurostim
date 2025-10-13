function [F0, dFF] = baseline_maximin(F, fs, sig_baseline, win_baseline,pv)
% BASELINE_MAXIMIN  Suite2p-style "maximin" baseline and ΔF/F
%   [F0, dFF] = baseline_maximin(F, fs, sig_baseline, win_baseline)
%
% Inputs
%   F            : T×N (double or single). Neuropil-corrected fluorescence.
%   fs           : sampling rate (Hz).
%   sig_baseline : Gaussian sigma in seconds (e.g., 10).
%   win_baseline : min/max window length in seconds (e.g., 60).
%
% 'casual' : use causal filters (defaults to true).
% 'pad_method' : method for padding data to reduce boundary effects.
%                Options: "reflect", "symmetric", "replicate", "none" (default: "reflect")
%
% Outputs
%   F0  : T×N baseline estimate (single precision)
%   dFF : T×N ΔF/F = (F - F0) ./ F0  (single precision)
%
% Notes
%   • In short: Gaussian smoothing (σ = sig_baseline*fs frames),
%     then moving MIN (window = win_baseline*fs frames),
%     then moving MAX with same window.
%   • Vectorized across ROIs; supports gpuArray if available.
%   • NaNs in F are propagated through the Gaussian; movmin/movmax ignore NaNs at edges.
%   • Data is padded to reduce boundary effects, then trimmed back to original length.
arguments
    F (:,:) single % Work in single precision for speed
    fs (1,1) single
    sig_baseline (1,1) single
    win_baseline (1,1) single
    pv.causal (1,1) logical = true
    pv.pad_method (1,1) string = "reflect"  % "reflect", "symmetric", "replicate", "none"
end

[T, N] = size(F);

% Calculate padding requirements
sig_frames = max(1, sig_baseline * fs);           % σ in frames
win_frames = max(1, round(win_baseline * fs));
khalf = ceil(3 * sig_frames);

% Total padding needed (conservative estimate)
pad_length = khalf + win_frames;

% Apply padding to reduce boundary effects
if pv.pad_method ~= "none"
    switch pv.pad_method
        case "reflect"
            % Reflect data at boundaries (excludes edge point)
            pad_start = min(pad_length, T-1);
            pad_end = min(pad_length, T-1);
            if pad_start > 0
                F_start_pad = flipud(F(2:pad_start+1, :));
            else
                F_start_pad = single.empty(0, N);
            end
            if pad_end > 0
                F_end_pad = flipud(F(T-pad_end:T-1, :));
            else
                F_end_pad = single.empty(0, N);
            end
            F_padded = [F_start_pad; F; F_end_pad];

        case "symmetric"
            % Symmetric reflection (includes edge point)
            pad_start = min(pad_length, T);
            pad_end = min(pad_length, T);
            if pad_start > 0
                F_start_pad = flipud(F(1:pad_start, :));
            else
                F_start_pad = single.empty(0, N);
            end
            if pad_end > 0
                F_end_pad = flipud(F(T-pad_end+1:T, :));
            else
                F_end_pad = single.empty(0, N);
            end
            F_padded = [F_start_pad; F; F_end_pad];

        case "replicate"
            % Extend edge values
            F_padded = [repmat(F(1,:), pad_length, 1); F; repmat(F(T,:), pad_length, 1)];
    end
    % Calculate the start index for extracting original data later
    pad_start_samples = size(F_padded, 1) - size(F, 1) - pad_length + 1;
else
    F_padded = F;
    pad_start_samples = 1;
end

% ---------- 1) Gaussian smoothing ----------
% Build a 1D Gaussian kernel (length ~ 6σ + 1)
khalf = ceil(3 * sig_frames);
x = (-khalf:khalf);
g = exp(-0.5 * (x ./ sig_frames).^2);
g = g / sum(g);                                   % normalize
g = single(g);

% Convolve along time (dimension 1), all ROIs in parallel
% conv2 is typically faster than filter/smoothdata for big matrices
F_smooth = conv2(F_padded, g', 'same');

% ---------- 2) Moving minimum ----------
win_frames = max(1, round(win_baseline * fs));
% movmin/movmax are fast and GPU-enabled in recent MATLAB
if pv.causal
    F_min = movmin(F_smooth, [win_frames-1 0], 1);    % causal window is a bit faster
else
    F_min = movmin(F_smooth, win_frames, 1);
end

% ---------- 3) Moving maximum ----------
if pv.causal
    F0_padded = movmax(F_min, [win_frames-1 0], 1);     % same window as above
else
    F0_padded = movmax(F_min, win_frames, 1);     % same window as above
end

% Extract original time range from padded result
if pv.pad_method ~= "none"
    F0 = F0_padded(pad_start_samples:pad_start_samples + T - 1, :);
else
    F0 = F0_padded;
end

if nargout>1
    dFF = (F - F0) ./ F0;
end
end
