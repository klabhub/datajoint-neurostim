classdef SegmentedTimeSeries < matlab.mixin.Copyable

    % SegmentedTimeSeries class for managing and extracting epochs from time series data.
    % This class facilitates the segmentation of a time series based on onset times,
    % allowing for the extraction of epochs (segments of data) around these onsets.
    % It also handles edge cases where the desired epoch window extends beyond the
    % boundaries of the original time series.

    properties

        sampling_time (1,:) double % (1, n_sample) timepoints of the original signal
        signal (:,:) double % (n_channel, n_sample) original signal

        onsets (:,1) double % (1, n_trial) onset times of trials/epochs
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

            t = self.sampling_time;
            self.iStart = arrayfun(@(i) gen.absargmin(t - (self.onsets(i) + self.epoch_window(1))), 1:self.n_trial);
            self.iOnset = self.iStart + find(self.timepoints >= 0, 1, 'first') - 1;
            self.iEnd = self.iStart + self.epoch_length-1;

        end

        function ep = get.raw_epochs(self)

            nTrials = self.n_trial;
            nChannels = self.n_channel;
            nEpochSamples = self.epoch_length;

            signalN = self.signal;         
           
            isBefore = self.iStart<=0;
            n_nanpad_presignal=0;
            if any(isBefore)

                n_nanpad_presignal = self.iStart(find(~isBefore,1,'first')) - self.iStart(1) + 1;
                signalN = [nan(nChannels,n_nanpad_presignal), signalN];

            end

            isAfter = self.iEnd>self.n_sample;
            if any(isAfter)

                n_nanpad_postsignal = self.iEnd(end) - self.n_sample;
                signalN = [signalN, nan(nChannels,n_nanpad_postsignal)];
                
                
            end
            
            ep = nan([nTrials, nChannels, nEpochSamples]);

            for iEp = 1: nTrials
                
                for iCh = 1:nChannels

                    ep(iEp, iCh, :) = signalN(iCh, n_nanpad_presignal + (self.iStart(iEp):self.iEnd(iEp)));
                
                end               

            end

        end

        function ep = get.epochs(self)

            ep = self.raw_epochs;

            if self.isDetrend

                ep = self.detrend_(ep,3);

            end

            if self.isBaseline

                isBaselineT = gen.ifwithin(self.timepoints,self.baseline_window);
                ep = ep - mean(ep(:,:,isBaselineT),3);

            end

        end

        function ep = get_raw_epochs(self)

            ep = self.raw_epochs;

        end

        function l = get.epoch_length(self)

            l = length(self.timepoints);

        end

        function t = get.timepoints(self)
            % Returns the time points corresponding to the extracted epochs.
            % Assumes all epochs have the same length.

            if isempty(self.resolution) || isempty(self.epoch_window)
                t = []; % Return empty if resolution or epoch_window are not set.
                return;
            end

            t = self.epoch_window(1):self.resolution:self.epoch_window(2); %Create time vector.

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

