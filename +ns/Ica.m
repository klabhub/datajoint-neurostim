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
            v = (ns.C * proj(ns.CParm)) * (ns.IcaParm & 'session=0');
        end
    end



    methods (Access=public)

        function plot(tbl, pv)
            % topoplots for ICA
            %
            % To select components, use a find struct (see ns.Label/find)
            %  f.ltag= "iclabel" - search in Label with ltag iclabel
            % Search components with specific labels:
            %  f.value = ["Eye" "Muscle" "Other"]; 
            % or search components where the probability for certain
            % artifacts is above some threshold value:
            % f.value = {["Eye" "Muscle" "Other"],0.8}; 
            % plot(tbl,find = f)
            arguments
                tbl (1,1) ns.Ica
                pv.comp (1,:) double {mustBeInteger,mustBePositive} = 1:12
                pv.labels (1,:) string = ""  % Which ns.Label ltag to use for labeling components. Defaults to all
                pv.find  = struct.empty   % Optional struct to filter components; passed to ns.Label(Session).find()
                pv.tilesPerFigure (1,1) double = 24
                pv.etaDisplay (1,1) string {mustBeMember(pv.etaDisplay, ["button", "inline", "none"])} = "button"
            end
            tpl = fetch(tbl, '*');
            assert(~isempty(tpl), 'No rows in ns.Ica to plot.');

            for i = 1:numel(tpl)
                this = tpl(i);
                labelRelvar = ns.Label  & this;
                if pv.labels~=""
                    labelRelvar = labelRelvar & struct('ltag', cellstr(pv.labels)');
                end
                compsToPlot = pv.comp;
                if ~isempty(pv.find)
                    findLabelRelvar =  ns.Label  & this  & (ns.LabelParm & struct('ltag',pv.find.ltag));
                    if ~isfield(pv.find,'op')
                        pv.find.op = function_handle.empty;
                    end
                    foundComps = find(findLabelRelvar ,pv.find.value,op=pv.find.op);
                    if isempty(foundComps)
                        continue;
                    else
                        compsToPlot = cat(1,foundComps.components{:});
                    end
                end

                labels = fetch(labelRelvar* ns.LabelParm, 'q', 'extra', 'parms');
                chanlocs = [fetch(ns.CChannel & this, 'channelinfo').channelinfo];
                w        = ns.Ica.getWeights(ns.stripToPrimary(ns.Ica,this));
                expName  = sprintf('%s @ %s %s | %s/%s', this.subject, this.session_date, this.starttime, this.ctag, this.itag);
                ns.Ica.plotComponents(w.winverse, this.variance, chanlocs, labels, compsToPlot, expName, pv.tilesPerFigure, pv.etaDisplay);
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
                    'itag', key.itag, ...
                    'ctag', key.ctag);

                src = fetch1(ns.IcaSession & sessionKey, ...
                        'nrcomponents', 'chanlabels', 'variance');
                % Remap session channel labels to per-experiment indices
                expChanlocs = fetch(ns.CChannel & key, 'channelinfo');
                expLabels = {expChanlocs.channelinfo.labels};
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
                % ns.C stores full data lenght with NaN for artifact
                % windows if clean_rawdata has been used. Remove those
                % windows for ICA.
                EEG = ephys.eeglab.removeNan(EEG);
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
                    'itag', key.itag,'ctag',key.ctag);
                src = fetch1(ns.IcaSession & sessionKey, ...
                    'winverse', 'sphere', 'weights', 'chanlabels');
                src.channels = fetch1(ns.Ica & key, 'channels');
                % chanlabels are channel label strings; caller must remap to
                % per-experiment indices using ismember on EEG.chanlocs.
            else
                src = fetch(ns.Ica & key, 'winverse', 'sphere', 'weights', 'channels');
                src.chanlabels = [];  % not needed for per-exp ICA
            end
        end

        function plotComponents(winverse, variance, chanlocs, labels, compsToPlot, expName, tilesPerFigure, etaDisplay)
            % Render ICA component topoplots in tiled figures.
            % winverse      - [nChans x nComps] mixing matrix
            % variance      - [1 x nComps] variance explained per component
            % chanlocs      - EEGLAB chanlocs struct array
            % labels        - struct array from fetch(ns.Label*ns.LabelParm,...)
            %                 with fields q, extra, parms (may be empty)
            % compsToPlot   - component indices to plot
            % expName       - string used as figure/sgtitle name
            % tilesPerFigure - number of tiles per figure
            % etaDisplay    - "button" for separate ETA figure, "inline" for per-tile inset, "none" to suppress ETA traces
            warning('off'); %'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout'
            figCntr    = 0;
            cCntr      = 0;
            currentFig = [];
            figName    = "";
            figTiles   = {};   % {compIdx, etaEntries} for every component in the current figure

            for c = compsToPlot(:)'
                if mod(cCntr, tilesPerFigure) == 0
                    % Finish previous figure: add ETA button if any tile had ETA data.
                    if etaDisplay == "button" && ~isempty(figTiles) && ~isempty(currentFig)
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
                % Use evalc to call topoplot to avoid the warnings that
                % cannot be turned off (due to nexttile and scaling in
                % topoplot).
                [~, ~] = evalc("topoplot(winverse(:, c), chanlocs, 'verbose', 'off', 'electrodes', 'off', 'numcontour', 8)");
                str = "#" + string(c) + " Var:" + string(round(variance(c), 1)) + "%";
                etaEntries = {};
                for l = 1:numel(labels)
                    switch labels(l).parms.method
                        case "iclabel"
                            str = [str; labels(l).q{c} + ":" + string(round(100*max(labels(l).extra(c,:), [], 2))) + "%"]; %#ok<AGROW>
                        case "eta"
                            if ~isempty(labels(l).q)
                                str = [str; labels(l).parms.plugin + ":" + strjoin(string(labels(l).parms.events), "/") + " D= " + string(round(labels(l).q(c), 2))]; %#ok<AGROW>

                                if ~isempty(labels(l).extra)
                                    nSamps = size(labels(l).extra, 2);
                                    tMs    = linspace(labels(l).parms.window(1), labels(l).parms.window(2), nSamps);
                                    etaEntries{end+1} = struct( ...  %#ok<AGROW>
                                        'timesMs', tMs, ...
                                        'eta',     labels(l).extra(c, :), ...
                                        'label',   labels(l).parms.plugin + ":" + strjoin(string(labels(l).parms.events), "/"));
                                end
                            end
                        case "spearman"
                            if ~isempty(labels(l).q)
                                str = [str; labels(l).parms.ctag + ":" + string(labels(l).parms.channel) + " r= " + string(round(labels(l).q(c), 2))]; %#ok<AGROW>
                            end
                        case "eog"
                            if ~isempty(labels(l).q)
                                str = [str; "EOG : r= " + string(round(labels(l).q(c), 2))]; %#ok<AGROW>
                            end
                    end
                end
                title(str);
                if etaDisplay == "inline"
                    ns.Ica.addEtaTile(c, etaEntries);
                end
                figTiles{end+1} = struct('compIdx', c, 'etaEntries', {etaEntries}); %#ok<AGROW>
            end
            % Finish the last figure.
            if etaDisplay == "button" && ~isempty(figTiles) && ~isempty(currentFig)
                ns.Ica.addEtaButton(currentFig, figTiles, figName);
            end
            sgtitle(expName);
        end

        function addEtaTile(compIdx, etaEntries)
            nexttile
            if isempty(etaEntries)
                axis off;
                title("#" + string(compIdx) + " ETA");
                return;
            end

            hold on;
            for j = 1:numel(etaEntries)
                e = etaEntries{j};
                plot(e.timesMs, e.eta, 'LineWidth', 1, 'DisplayName', e.label);
            end
            xline(0, 'k:');
            hold off;
            grid on;
            box off;
            xlabel('Time (ms)');
            title("#" + string(compIdx) + " ETA");
            if numel(etaEntries) > 1
                legend('Location', 'best');
            end
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

        function [signal,info] = clean(signal,cTpl,pv)
            arguments
                signal (:,:) double
                cTpl (1,1) struct
                pv.itag (1,1) string
                pv.ltag (1,1) string
                pv.find (1,1) struct= struct.empty
            end

            % ICA based preprocessing
            icaRelVar = ns.Ica & (ns.C & cTpl) & struct('itag',pv.itag);
            labelRelVar  = ns.Label & icaRelVar & struct('ltag',pv.ltag);
            if ~exists(icaRelVar)
                % Try a session level ICA
                icaRelVar = ns.IcaSession & (ns.Session & (ns.C &  cTpl)) & struct('itag',pv.itag);
                labelRelVar  = ns.LabelSession & icaRelVar & struct('ltag',pv.ltag);
            assert(exists(icaRelVar),"No ica with itag %s found for %s@%sT%s",pv.itag,cTpl.subject,cTpl.session_date,cTpl.starttime);
            W = ns.Ica.getWeights(fetch(icaRelVar));
            assert(exists(labelRelVar),"No Label found for %s in %s",pv.ltag,pv.itag);
            findLabelRelvar =  ns.Label  & icaRelVar  & (ns.LabelParm & struct('ltag',pv.ltag));
            if ~isfield(pv.find,'op')
                pv.find.op = function_handle.empty;
            end
            T = find(findLabelRelvar ,pv.find.value,op=pv.find.op);              
            compsToRemove = [T.components{:}];
            % Reconstruct the signal from these components.
            Xica = signal(:, W.channels);
            A = W.weights * W.sphere * Xica';
            Xclean = Xica' - W.winverse(:, compsToRemove) * A(compsToRemove, :);
            signal(:, W.channels) = Xclean';
            varExplained = fetch1(icaRelVar,'variance');
            info.nrComponents= numel(compsToRemove);
            info.variance = sum(varExplained(compsToRemove));
            fprintf('Removed %d components (%.0f%% variance) based on ICA %s\n',info.nrComponents,info.variance,pv.itag);    
        end
    end

end
