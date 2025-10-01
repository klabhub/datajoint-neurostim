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
arguments
    F (:,:) single % Work in single precision for speed 
    fs (1,1) single 
    sig_baseline (1,1) single
    win_baseline (1,1) single
    pv.causal (1,1) logical = true
end 
   
    % ---------- 1) Gaussian smoothing ----------
    sig_frames = max(1, sig_baseline * fs);           % σ in frames
    % Build a 1D Gaussian kernel (length ~ 6σ + 1)
    khalf = ceil(3 * sig_frames);
    x = (-khalf:khalf);
    g = exp(-0.5 * (x ./ sig_frames).^2);
    g = g / sum(g);                                   % normalize
    g = single(g);

    % Convolve along time (dimension 1), all ROIs in parallel
    % conv2 is typically faster than filter/smoothdata for big matrices
    F_smooth = conv2(F, g', 'same');

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
        F0 = movmax(F_min, [win_frames-1 0], 1);     % same window as above
    else
        F0 = movmax(F_min, win_frames, 1);     % same window as above
    end

    if nargout>1 
        dFF = (F - F0) ./ F0;
    end
end
