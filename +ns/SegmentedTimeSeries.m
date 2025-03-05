classdef SegmentedTimeSeries < matlab.mixin.Copyable

    % SegmentedTimeSeries class for managing and extracting epochs from time series data.
    % This class facilitates the segmentation of a time series based on onset times,
    % allowing for the extraction of epochs (segments of data) around these onsets.
    % It also handles edge cases where the desired epoch window extends beyond the 
    % boundaries of the original time series.

    properties

        sampling_time (1,:) double {mustBeRow} % (1, n_sample) timepoints of the original signal
        signal (:,:) double {mustBeMatrix} % (n_channel, n_sample) original signal

        onsets (:,1) double {mustBeColumn} % (1, n_trial) onset times of trials/epochs
        epoch_window (1,2) double % Time window [t_before, t_after] relative to each onset.

        iOnset (:,1) int64 % Indices of the onset times within the timepoints array.
        iStart (:,1) int64 % Start indices of the epochs.
        iEnd (:,1) int64  % End indices of the epochs.

        trial_labels % Labels associated with each trial (optional).
        baseline_window = []

    end

    properties (Access = private)

        isBaseline = false
        isDetrend = false

    end

    properties (Access = private, Dependent)

        raw_epochs

    end


    properties (Dependent)

        n_trial % Number of trials/segments.
        n_sample % Number of samples in the original time series.
        n_channel % Number of channels in the time series data.
        resolution % Time resolution of the time series (assumed regular spacing).
        timepoints % Time points corresponding to the extracted epochs.
        epoch_length % Length of an epoch
        isMissingData % Flag indicating if missing data was encountered during epoch extraction.
        epochs % Extracted epochs (n_trial x n_epoch_samples x n_channel).

    end

    methods

        function self = SegmentedTimeSeries(t, t_onsets, epoch_win, signal, labels)

            self.signal = signal;
            self.sampling_time = t;
            self.onsets = t_onsets;
            self.epoch_window = epoch_win;
            if nargin > 4, self.trial_labels = labels; end

        end

        % Get functions for dependent properties
        function n = get.n_trial(self)
            n = length(self.onsets);
        end

        function n = get.n_sample(self)
            n = length(self.sampling_time);
        end

        function n = get.n_channel(self)
            [n, ~] = size(self.signal);  % Get the size; n will be the number of rows (channels)
        end

        function res = get.resolution(self)
            res = self.sampling_time(2) - self.sampling_time(1); % Assumes regular spacing
        end


        function self = make(self)

            [iOnsets, iStarts, iEnds] = self.segment_indices(self.sampling_time, self.onsets, self.epoch_window);
            self.iOnset = iOnsets;
            self.iStart = iStarts;
            self.iEnd = iEnds;

        end

        function ep = get.raw_epochs(self)
            
            nTrials = self.n_trial;
            nChannels = self.n_channel;
            nEpochSamples = max(self.iEnd - self.iStart + 1);

            % Create all indices at once (nTrials x nEpochSamples)
            all_indices = repmat(self.iStart, 1, nEpochSamples) + repmat((0:nEpochSamples-1), nTrials, 1);

            % Bounds check for all indices (nTrials x nEpochSamples)
            valid_indices = (all_indices >= 1) & (all_indices <= self.n_sample);
            self.signal(:, all_indices(~valid_indices)) = NaN;
            ep = permute( ...
                reshape(...
                self.signal(:,all_indices), ...
                nChannels, nTrials, nEpochSamples),...
                [2,1,3]);
            self.signal(:, all_indices(~valid_indices)) = [];

        end

        function ep = get.epochs(self)

            ep = self.raw_epochs;

            if self.isDetrend

                ep = self.detrend_(ep,3);

            end

            if self.isBaseline
                
                ep = ep - mean(self.extract_timepoints(ep, self.timepoints, self.baseline_window), 3, 'omitnan');

            end

        end

        function ep = get_raw_epochs(self)

            ep = self.raw_epochs;

        end

        function l = get.epoch_length(self)

            l = size(self.raw_epochs, 3);

        end
        
        function t = get.timepoints(self)
            % Returns the time points corresponding to the extracted epochs.
            % Assumes all epochs have the same length.

            if isempty(self.resolution) || isempty(self.epoch_window)
                t = []; % Return empty if resolution or epoch_window are not set.
                return;
            end

            t_start = self.epoch_window(1);  % Time relative to onset for the start of the epoch.

            t = t_start + (0:self.epoch_length-1) * self.resolution; %Create time vector.

        end

        function i = get.isMissingData(self)

            i = any(isnan(self.raw_epochs),[2,3]);
            
        end


        function baseline(self, t_win)

            self.isBaseline = true;
            self.baseline_window = t_win;

        end

        function detrend(self)

            self.isDetrend = true;
        end

    end

    methods (Static, Access = protected)

        function [iOnsets, iStarts, iEnds] = segment_indices(t, t_onsets, t_win)
            
            arguments
                t (1,:) double %timeseries
                t_onsets (:,1) double % onset time of each segment, t0
                t_win (1, 2) double % time window [t_before, t_after] relative to t0
            end

            nOnsets = length(t_onsets);
            dt = t(2) - t(1); % Calculate the time step (assuming regular spacing)

            iOnsets = zeros(nOnsets,1);
            iStarts = zeros(nOnsets,1);
            iEnds = zeros(nOnsets,1);

            for i = 1:nOnsets
                t0 = t_onsets(i);
                [~, iOnset] = min(abs(t - t0));
                iOnsets(i) = iOnset;

                % t_start = t0 + t_win(1);
                % [~, iStart] = min(abs(t - t_start));
                % 
                % t_end = t0 + t_win(2);
                % [~, iEnd] = min(abs(t - t_end));

                iStart = round(iOnset + t_win(1)/dt); % Calculate start index directly
                iEnd = round(iOnset + (t_win(2)/dt))-1;  % Calculate end index directly

                % Handle edge cases
                if iStart < 1
                    warning("Onset %d: Start of window outside timeseries, adjusting", i);
                    iStart = 1;
                end
                if iEnd > length(t)
                    warning("Onset %d: End of window outside timeseries, adjusting", i);
                    iEnd = length(t);
                end

                iStarts(i) = iStart;
                iEnds(i) = iEnd;
            end
        end

        function s = extract_timepoints(ep, t, t_win)

            % Extracts data at given time points
            isT = (t_win(1) <= t) & (t_win(2) >= t);
            s = ep(:, :, isT);

        end

        function A = detrend_(A, dim)

            % A: Input array
    
            % Get the number of dimensions
            n_dims = ndims(A);
            if n_dims == 2, A = detrend(A); return; end
        
            % Create the permutation order
            perm_order = 1:n_dims;
            perm_order([dim, n_dims]) = perm_order([n_dims, dim]);
        
            % Permute the array
            A = permute(A, perm_order);
        
            % Detrend along the last dimension
            A = detrend(A);
        
            % Permute back to the original order
            A = permute(A, perm_order);

        end
    end

end

