%{
# Preprocessed epoch data, the actual epochs are stored in ns.EpochChannel

-> ns.C
-> ns.EpochParm
trial : int
---
onset : float
flag : varchar(32)
%}

classdef Epoch < dj.Computed & dj.DJInstance

    properties

        data = []
        frequencies = []

    end

    properties (Dependent)

        C % ns.C table
        Experiment % ns.Experiment table
        EpochParm % ns.EpochParm table
        EpochChannel

        channels
        timepoints   
        conditions

        keySource

    end

    methods

        function v = get.keySource(~)

            v = ns.C * ns.Experiment * proj(ns.EpochParm,'dimension');
        
        end
    end

    methods (Access = protected)

        function makeTuples(eTbl, key)

            expTbl = ns.Experiment & key;           

            epTbl = ns.EpochParm & key;
            cTbl = ns.C & key;

            abs_onsets = expTbl.first_frame_onsets;
            rel_onsets = get(expTbl, epTbl{'plugin'}, 'prm', 'startTime', 'what', 'trialtime');
                        
            if isempty(abs_onsets), return; end
            

            % trials that were shown

            isVld = ~isinf(rel_onsets) & ~isnan(rel_onsets);
            trial_no = get(expTbl, epTbl{'plugin'}, 'prm', 'startTime', 'what', 'trial');
            trial_no = trial_no(isVld);
            onsets = abs_onsets(trial_no) + rel_onsets(trial_no);
            
            t = cTbl.timepoints;
            dt = cTbl.dt;
            ep_win = epTbl{'epoch_win'};
            ep_timepoints = ep_win(1):dt:ep_win(2);

            if isempty(trial_no)
                fprintf("No epochs were detected for:\n");
                disp(key);
                return
            else
                fprintf("\t%d valid trials were found.", length(trial_no));
            end

            t_export = gen.Timer().start("Exporting data from ns.CChannel...\n");

            signal = fetch(ns.CChannel & cTbl, 'signal');
            signal = horzcat(signal(:).signal)';

            t_export.stop("\tExporting is complete.");
            t_export.report();

            t_sgm = gen.Timer().start("Now segmenting...\n");
            
            ep = segment();
            

            t_sgm.stop("\t\tSegmentation complete.");
            t_sgm.report();

            pv = epTbl{'pv'};
            
            if pv.detrend

                ep = eTbl.detrend(ep);

            end

            if pv.baseline

                ep = eTbl.baseline(ep, gen.ifwithin(ep_timepoints, pv.baseline_win));

            end

            if pv.rereference

                ep(:, cTbl.artifacts.all,:) = NaN; % exclude noisy channels
                ep = eTbl.rereference(ep, cTbl.channel_coordinates, pv.ref_opts{:});

            end

            n_art_fun = length(pv.artifact_parm);
            flags = cell(1,n_art_fun);
            excludeN = [];
            epoch_flags = repmat("", length(isVld), 1);

            for iFun = 1:n_art_fun
                
                fun_str = pv.artifact_parm(iFun).fun;

                t_artN = gen.Timer().start("Finding artifacts through: %s.", fun_str);
                if strcmp(fun_str, 'detect_outliers')

                    %make sure handle refers to the nested function

                    funN = @detect_outliers;

                else

                    funN = str2func(fun_str);
                end
                
                argN = pv.artifact_parm(iFun).args;
                if ~iscell(argN)
                    % if the input arguments are not stored in a cell array
                    % transform it to cell to make compatible for varargin
                    argN = {argN};
                end
                
                flags{iFun} = funN(key, 'exclude', excludeN, argN{:});
                % exclude bad epochs from the prior step
                excludeN = [excludeN; flags{iFun}.all];
                % update epoch flags
                epoch_flags = flags{iFun}.flag(epoch_flags);
                t_artN.stop("\tArtifact detection '%s' complete.", fun_str);
                t_artN.report();

            end

            %% Submission
            t_sub = gen.Timer().start("\tNow submitting to the server\n");

            % Create epoch tuple (n_epochs,1)
            epoch_tpl = mergestruct(key, struct( ...
                trial = num2cell(trial_no),...
                onset = num2cell(rel_onsets(isVld)),...
                flag = num2cell(epoch_flags(isVld))...
                ));
            chunkedInsert(ns.Epoch, epoch_tpl);

            [n_epochs, n_channels, n_timepoints] = size(ep);
            
            % reshape epoch to flatten the first two dimensions, channels
            % unravel first, then epochs
            ep = reshape(permute(ep, [2, 1, 3]), n_epochs*n_channels, n_timepoints);
            
            nonprimKey = setdiff(fieldnames(epoch_tpl), eTbl.primaryKey);
            epoch_tpl = mergestruct( ...
                repelem(rmfield(epoch_tpl, nonprimKey), n_channels),...
                struct(channel = repmat(gen.make_column(num2cell(cTbl.channels)), n_epochs, 1),...
                signal = num2cell(ep, 2)...
            ));
            chunkedInsert(ns.EpochChannel, epoch_tpl)

            t_sub.stop("\t\tSubmission is complete.\n");
            t_sub.report();

            function ep = segment()
                
                % Get indices for each epoch start, onset, and end
                
                nTrials = length(onsets);
                nChannels = size(signal,1);
                nEpochSamples = length(ep_timepoints);
                nSample = size(t, 2);
                iStart = arrayfun(@(i) gen.absargmin(t - (onsets(i) + ep_win(1))), 1:nTrials);
                % iOnset = iStart + find(ep_timepoints >= 0, 1, 'first') - 1;
                iEnd = iStart + nEpochSamples - 1;
               
                signalN = signal;

                isBefore = iStart <= 0;
                n_nanpad_presignal = 0;
                if any(isBefore)

                    n_nanpad_presignal = iStart(find(~isBefore,1,'first')) - iStart(1);
                    signalN = [nan(nChannels, n_nanpad_presignal), signalN];

                end

                isAfter = iEnd > nSample;
                if any(isAfter)

                    n_nanpad_postsignal = iEnd(end) - nSample;
                    signalN = [signalN, nan(nChannels, n_nanpad_postsignal)];

                end

                ep = nan([nTrials, nChannels, nEpochSamples]);

                for iEp = 1: nTrials

                    for iCh = 1:nChannels

                        ep(iEp, iCh, :) = signalN(iCh, n_nanpad_presignal + (iStart(iEp):iEnd(iEp)));

                    end

                end


            end

            function flags = detect_outliers(varargin)
                
                flags = ns.detect_outlier_epochs(ep, cTbl.srate, 'epoch_no', trial_no, varargin{2:end});

            end            
            
        end
        
    end
    
    methods

        function varargout = transform(eTbl, varargin)
            % Delete this in future, now all is done within EpochChannel
                                
                cTbl = eTbl.C;
                epTbl = eTbl.EpochParm;
                etags = unique(epTbl{'etag'},'rows');
                ctags = unique(cTbl{'ctag'},'rows');
                files = unique(cTbl{'filename'}, 'rows');


                if size(etags,1)>1 || size(ctags,1)>1 || size(files,1) > 1
                    error("Transformed data contains data from different files and/or preprocessing pipelines.");
                end

                if isempty(eTbl.data)

                    ep = eTbl.load();

                else

                    ep = eTbl.data;

                end
                n_arg = nargin - 1;
                ii = 1;

                while ii <= n_arg

                    argN = varargin{ii};
                    ii = ii + 1;

                    switch argN
                        case 'include'

                            argN = varargin{ii};
                            ep = eTbl.include(ep, argN{:});
                            
                            ii = ii + 1;
                            if ii > n_arg

                                varargout{1} = ep;

                            end
                            
                        case 'exclude'
                        case 'baseline'
                        case 'detrend'
                        case 'fft'

                            [ep, ph, fq] = eTbl.do_fft(ep, cTbl.srate, 3);
                            eTbl.frequencies = fq;
                            if ii > n_arg

                                varargout{1} = ep;
                                varargout{2} = ph;
                                varargout{3} = fq;

                            end
                            
                        case 'snr'    

                            argN = varargin{ii}; %options
                            ii = ii + 1;

                            fq_res = eTbl.frequencies(2)-eTbl.frequencies(1);
                            fq_bin_skip = argN{1};
                            fq_bin = argN{2};
                            
                            n_bins = round(fq_bin/fq_res);
                            n_bins_skip = round(fq_bin_skip/fq_res);
                            mid_bin = n_bins + 1;
                            kernel = true(1, 2*n_bins + 1);
                            kernel(mid_bin-n_bins_skip : mid_bin+n_bins_skip) = false;

                            if numel(argN) > 2
                                method = argN{3};
                            else
                                method = @mean;
                            end

                            noise = eTbl.calculate_noise_(ep, kernel, 3, method);
                            ep = ep./noise;

                            if ii > n_arg
                                varargout{1} = ep;
                                varargout{2} = eTbl.frequencies;
                            end

                    end

                end

                eTbl.data = ep;
                
        end

        function ep = load(eTbl, qry)

            % loads data
            if nargin < 2 || isempty(qry)

                epTbl = eTbl.EpochChannel;

            else

                epTbl = eTbl.EpochChannel & qry;          

            end
            %% Fetch data                
            t_fetch = gen.Timer().start('Fetching epochs...');

            
            ep = fetch(epTbl , 'signal');
            ep = cat(1, ep(:).signal);

            % get dimension size of the final ep matrix
            n_ep = count(eTbl);
            n_ch = count(ns.CChannel & eTbl.C);
            n_t = size(ep,2);

            % reshape ep to n_ch by n_ep by n_t
            ep = reshape(ep, n_ch, n_ep, n_t);
            %permute to get the desired shape n_ep n_ch, n_t
            ep = permute(ep, [2, 1, 3]);

            eTbl.data = ep;

            t_fetch.stop("complete.\n");
            t_fetch.report();

        end

        function ep = include(eTbl, ep, pv)

            arguments

                eTbl ns.Epoch
                ep double
                pv.channels = []
                pv.time_window = []
                
                pv.trials = []
                pv.conditions = []
            
            end

            if ~isempty(pv.channels)

                isChIn = ismember(eTbl.channels, pv.channels);
                ep = ep(:,isChIn,:);

            end

            if ~isempty(pv.time_window)

                isTIn = gen.ifwithin(eTbl.timepoints, pv.time_window);
                ep = ep(:,:,isTIn);

            end

            if ~isempty(pv.trials)

                trl = fetch(eTbl,'trial');
                trl = [trl(:).trial];
                isTrlIn = ismember(trl, pv.trials);

                ep = ep(isTrlIn, :, :);

            end

            if ~isempty(pv.conditions)

                isTrlIn = ismember(eTbl.conditions, pv.conditions);
                ep = ep(isTrlIn, :, :);

            end
            


        end

        function c = get.C(eTbl)

            c = ns.C & eTbl;

        end

        function exp = get.Experiment(eTbl)

            exp = ns.Experiment & eTbl;

        end

        function e = get.EpochParm(eTbl)

            e = ns.EpochParm & eTbl;

        end

        function e = get.EpochChannel(eTbl)

            e = ns.EpochChannel & eTbl;
            
        end


        function ch = get.channels(eTbl)

            ch = eTbl.C.channels;
            
        end

        function cond = get.conditions(eTbl)
            
            dTbl = ns.DimensionCondition & eTbl;
            dTbl = fetch(dTbl,'trials');
            trials = fetch(eTbl,'trial');
            trials = [trials(:).trial]';

            cond = repelem("",length(trials),1);
            for ii = 1:length(dTbl)

              
                cond(ismember(trials,dTbl(ii).trials)) = dTbl(ii).name;
            
            end

        end

        function t = get.timepoints(eTbl)

            cTbl = ns.C & eTbl;

            epTbl = ns.EpochParm & eTbl;
            ep_win = epTbl{'epoch_win'};
            dt = cTbl.dt;
            t = ep_win(1):dt:ep_win(2);

        end
      
                
    end  

    methods (Static)

        function ep = detrend(ep, varargin)

            n_dim = ndims(ep);
            ep = permute(ep, circshift(1:n_dim,1));
            ep = detrend(ep,1,'omitnan', varargin{:});
            ep = permute(ep, circshift(1:n_dim,-1));

        end

        function ep = baseline(ep, isEpWin)

            base = mean(ep(:,:,isEpWin), ndims(ep), 'omitnan');
            ep = ep - base;

        end

        function ep = rereference(ep, chanLocs, varargin)

            varargin = [varargin, 'channel_locations', chanLocs];

            ep = ns.rereference(ep, varargin{:});

        end

        function [amplitude, phase, frequencies] = do_fft(data, fs, n)
            % fftSliceReal - Computes FFT amplitude and phase for each slice along the n'th dimension.
            %                Only includes real frequencies.
            %
            % Inputs:
            %   data: Input multi-dimensional data.
            %   n: Dimension along which to slice and compute FFT.
            %   fs: Sampling frequency (optional, default is 1).
            %
            % Outputs:
            %   amplitude: Amplitude of the FFT.
            %   phase: Phase of the FFT.
            %   frequencies: Corresponding real frequencies.

            n_dim = ndims(data);

            % Compute FFT for each slice
            fftResult = fft(data, [], n_dim);

            idx = repmat({':'},1,n_dim);
            % Calculate real frequencies
            N = size(data, n);
            if mod(N, 2) == 0
                frequencies = (0:N/2) * fs / N;
                idx{n_dim} = 1:N/2+1;
            else
                frequencies = (0:(N-1)/2) * fs / N;
                idx{n_dim} = 1:(N+1)/2;
            end

            amplitude = 2*abs(fftResult(idx{:})/sqrt(N));
            phase = angle(fftResult(idx{:}));


        end

        function noise_arr = calculate_noise_(arr, filter_kernel, n_dim, method)

            if nargin < 3 || isempty(n_dim), n_dim = ndims(arr); end
            if nargin < 4; method = @mean; end

            idx_selector = repmat({':'}, 1, ndims(arr));
            idx_assigner = repmat({':'}, 1, ndims(arr));
            
            n_points = size(arr, n_dim);
            noise_arr = zeros(size(arr));
            n_kernel = numel(filter_kernel);

            for i = 1:n_points
                idx_assigner{n_dim} = i;

                if i <= n_kernel/2 || i > n_points - n_kernel/2

                    noise_arr(idx_assigner{:}) = NaN;

                else

                    idxN = false(1, size(arr,n_dim));
                    idxN(floor(i+1-n_kernel/2):floor(i+n_kernel/2)) = filter_kernel;

                    idx_selector{n_dim} = idxN;

                    noise_arr(idx_assigner{:}) = method(arr(idx_selector{:}), n_dim);

                end

            end

        end

    end

end

