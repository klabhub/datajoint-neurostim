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

    properties (Dependent)

        C % ns.C table
        Experiment % ns.Experiment table

        channels
        timepoints
        frequencies    
        conditions

        keySource

    end

    properties (Access = protected)


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
            dims = fetch(ns.DimensionCondition & key, '*');

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

            t_export = gen.Timer().start("Exporting data from ns.CChannel...\n");

            signal = fetch(ns.CChannel & cTbl, 'signal');
            signal = horzcat(signal(:).signal)';
            ch_list = cTbl.channels;

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

        function ep = retrieve(eTbl, varargin)

            ep = ns.RetrievedEpochs(eTbl, varargin{:});

        end


        

    end

    methods % GET Methods

        function c = get.C(eTbl)

            c = ns.C & eTbl;

        end

        function exp = get.Experiment(eTbl)

            exp = ns.Experiment & eTbl;

        end


        function ch = get.channels(eTbl)

            ch = eTbl.C.channels;
            
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

    end

end

