%{
# Stores Independent Component Analysis results for a row in ns.C
-> ns.C
-> ns.IcaParm
---
nrcomponents :  int       # Number of components
winverse = NULL : longblob  # Inverse ICA weights (null for session-ICA; stored in ns.IcaSession)
sphere   = NULL : longblob  # Sphering matrix (null for session-ICA; stored in ns.IcaSession)
weights  = NULL : longblob  # ICA weights (null for session-ICA; stored in ns.IcaSession)
channels : longblob  # Channels used for ICA
variance : blob      # Variance explained per component
%}
%
% Currently implemented with EEGLAB
%
% See also ns.IcaParm, ephys.eeglab.preprocess, ns.Label, ns.LabelParm
%
% BK - May 2026
classdef Ica < dj.Computed & dj.DJInstance
    properties (Dependent)
        keySource
    end

    methods
        function v = get.keySource(~)
            % Per-experiment ICA: ns.C paired with non-session IcaParms.
            v = (ns.C * ns.CParm) * (ns.IcaParm & 'session=0');
        end
    end

    

    methods (Access=public)

        function plot(tbl, pv)
            arguments
                tbl (1,1) ns.Ica
                pv.comp (1,:) double {mustBeInteger,mustBePositive} = 1:12
                pv.labels (1,:) string = "iclabel"  % Which ns.Label ltag to use for labeling components in the plot
                pv.find  = string.empty % Optional string to filter components; e.g. find="blink"
                pv.tilesPerFigure (1,1) double = 24
            end
            tpl = fetch(tbl, '*');
            assert(~isempty(tpl), 'No rows in ns.Ica to plot.');

            for i = 1:numel(tpl)
                this = tpl(i);
                compsToPlot = pv.comp;
                if ~isempty(pv.find)
                    foundComps = find(ns.Label & this, pv.find);
                    if isempty(foundComps)
                        warning('No components found that match the find instruction. Plotting all components instead.');
                    else
                        compsToPlot = [foundComps.components{:}];
                    end
                end
                labels   = fetch(ns.Label*ns.LabelParm & this & struct('ltag', cellstr(pv.labels)'), 'q', 'extra', 'parms');
                chanlocs = [fetch(ns.CChannel & this, 'channelinfo').channelinfo];
                w        = ns.Ica.getWeights(this);
                expName  = sprintf('%s @ %s %s | %s/%s', this.subject, this.session_date, this.starttime, this.ctag, this.itag);
                ns.Ica.plotComponents(w.winverse, this.variance, chanlocs, labels, compsToPlot, expName, pv.tilesPerFigure);
            end
        end
    end

    methods (Access=protected)
        function makeTuples(tbl, key)
            isSession = fetch1(ns.IcaParm & key, 'session');

            if isSession
                % Weights are stored in ns.IcaSession — do not duplicate them here.
                % Only store lightweight per-experiment metadata.
                sessionKey = struct('subject', key.subject, ...
                                   'session_date', key.session_date, ...
                                   'itag', key.itag);
                src = fetch1(ns.IcaSession & sessionKey, ...
                    'nrcomponents', 'chanlabels', 'variance');
                % Remap session channel labels to per-experiment indices
                expChanlocs = fetch(ns.CChannel & key, 'channelinfo');
                %expLabels = {[expChanlocs.channelinfo].labels};
                [~, chanIdx] = ismember(src.chanlabels, expLabels);
                chanIdx = chanIdx(chanIdx > 0);
                tpl = key;
                tpl.nrcomponents = src.nrcomponents;
                tpl.channels     = chanIdx;   % per-experiment indices for this row
                tpl.variance     = src.variance;
                % winverse, sphere, weights left NULL — use ns.Ica.getWeights(key)
            else
                % Per-experiment ICA
                icaParms = fetch1(ns.IcaParm & key, 'parms');
                EEG = ephys.eeglab.dataset(key, data=key.ctag, itag="");
                if isfield(icaParms, 'filt')
                    EEG = pop_eegfiltnew(EEG, 'hicutoff', icaParms.filt.hicutoff, ...
                                              'locutoff',  icaParms.filt.locutoff);
                    icaParms = rmfield(icaParms, 'filt');
                end
                parms.eeglab.ica = icaParms;
                EEG = ephys.eeglab.preprocess(EEG, parms);

                varExplained = ephys.eeglab.icaVarianceExplained(EEG);

                tpl = key;
                tpl.nrcomponents = size(EEG.icaweights, 1);
                tpl.winverse     = EEG.icawinv;
                tpl.sphere       = EEG.icasphere;
                tpl.weights      = EEG.icaweights;
                tpl.channels     = EEG.icachansind;
                tpl.variance     = varExplained;
            end
            insert(tbl, tpl);
        end
    end

    methods (Static)
        function src = getWeights(key)
            % Fetch ICA weight matrices for a given ns.Ica key.
            % For session-scoped ICA (IcaParm.session=true), the weights are
            % stored once in ns.IcaSession and not duplicated in ns.Ica.
            % For per-experiment ICA they are in ns.Ica itself.
            % Returns a struct with fields: winverse, sphere, weights, channels.
            isSession = fetch1(ns.IcaParm & key, 'session');
            if isSession
                sessionKey = struct('subject', key.subject, ...
                                   'session_date', key.session_date, ...
                                   'itag', key.itag);
                src = fetch(ns.IcaSession & sessionKey, ...
                    'winverse', 'sphere', 'weights', 'chanlabels');
                % chanlabels are channel label strings; caller must remap to
                % per-experiment indices using ismember on EEG.chanlocs.
            else
                src = fetch(ns.Ica & key, 'winverse', 'sphere', 'weights', 'channels');
                src.chanlabels = [];  % not needed for per-exp ICA
            end
        end

        function plotComponents(winverse, variance, chanlocs, labels, compsToPlot, expName, tilesPerFigure)
            % Render ICA component topoplots in tiled figures.
            % winverse      - [nChans x nComps] mixing matrix
            % variance      - [1 x nComps] variance explained per component
            % chanlocs      - EEGLAB chanlocs struct array
            % labels        - struct array from fetch(ns.Label*ns.LabelParm,...)
            %                 with fields q, extra, parms (may be empty)
            % compsToPlot   - component indices to plot
            % expName       - string used as figure/sgtitle name
            % tilesPerFigure - number of tiles per figure
            warning('off'); %'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout'
            figCntr    = 0;
            cCntr      = 0;
            currentFig = [];
            figName    = "";
            figTiles   = {};   % {compIdx, etaEntries} for every tile in the current figure

            for c = compsToPlot(:)'
                if mod(cCntr, tilesPerFigure) == 0
                    % Finish previous figure: add ETA button if any tile had ETA data.
                    if ~isempty(figTiles) && ~isempty(currentFig)
                        ns.Ica.addEtaButton(currentFig, figTiles, figName);
                    end
                    figCntr    = figCntr + 1;
                    figName    = expName + "-" + string(figCntr);
                    currentFig = figByName(figName);
                    clf;
                    tiledlayout('flow');
                    figTiles = {};
                end
                cCntr = cCntr + 1;
                nexttile
                topoplot(winverse(:, c), chanlocs, 'verbose', 'off', 'electrodes', 'off', 'numcontour', 8);
                str = "#" + string(c) + " Var:" + string(round(variance(c), 1)) + "%";
                etaEntries = {};
                for l = 1:numel(labels)
                    switch labels(l).parms.method
                        case "iclabel"
                            str = [str; labels(l).q{c} + ":" + string(round(100*max(labels(l).extra(c,:), [], 2))) + "%"]; %#ok<AGROW>
                        case "eta"
                            str = [str; labels(l).parms.plugin + ":" + strjoin(string(labels(l).parms.events), "/") + " z= " + string(round(labels(l).q(c), 2))]; %#ok<AGROW>
                            if ~isempty(labels(l).extra)
                                nSamps = size(labels(l).extra, 2);
                                tMs    = linspace(-labels(l).parms.window, labels(l).parms.window, nSamps);
                                etaEntries{end+1} = struct( ...  %#ok<AGROW>
                                    'timesMs', tMs, ...
                                    'eta',     labels(l).extra(c, :), ...
                                    'label',   labels(l).parms.plugin + ":" + strjoin(string(labels(l).parms.events), "/"));
                            end
                        case "spearman"
                            str = [str; labels(l).parms.ctag + ":" + string(labels(l).parms.channel) + " r= " + string(round(labels(l).q(c), 2))]; %#ok<AGROW>
                    end
                end
                title(str);
                figTiles{end+1} = struct('compIdx', c, 'etaEntries', {etaEntries}); %#ok<AGROW>
            end
            % Finish the last figure.
            if ~isempty(figTiles) && ~isempty(currentFig)
                ns.Ica.addEtaButton(currentFig, figTiles, figName);
            end
            sgtitle(expName);
        end

        function addEtaButton(fig, figTiles, figName)
            % Add an "ETA" button in the top-right corner of fig if any tile has ETA data.
            hasEta = any(cellfun(@(t) ~isempty(t.etaEntries), figTiles));
            if ~hasEta, return; end
            uicontrol(fig, ...
                'Style',    'pushbutton', ...
                'String',   'ETA', ...
                'Units',    'normalized', ...
                'Position', [0.88 0.955 0.10 0.04], ...
                'Callback', @(~,~) ns.Ica.showEtaFigure(figTiles, figName));
        end

        function showEtaFigure(figTiles, figName)
            % Open (or refresh) a figure showing ETAs for every tile.
            % Tiles with no ETA are left blank to preserve positional correspondence.
            etaFig = figByName(figName + " eta");
            clf(etaFig);
            tiledlayout(etaFig, 'flow');
            for i = 1:numel(figTiles)
                d = figTiles{i};
                nexttile;
                if isempty(d.etaEntries)
                    axis off;
                    title("#" + string(d.compIdx));
                else
                    hold on;
                    for j = 1:numel(d.etaEntries)
                        e = d.etaEntries{j};
                        plot(e.timesMs, e.eta, 'DisplayName', e.label);
                    end
                    xline(0, 'k--', 'HandleVisibility', 'off');
                    xlabel('Time (ms)');
                    if numel(d.etaEntries) > 1, legend; end
                    title("#" + string(d.compIdx));
                end
            end
            sgtitle(figName + " - ETA");
        end
    end
 
end
