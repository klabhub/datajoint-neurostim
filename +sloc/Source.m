%{
# Source localized ssveps
-> ns.Epoch
-> sloc.Leadfield
-> sloc.Headmodel
-> sloc.SourceParm
label : varchar(32) # Some label to identify different conditions in an epoch
---
power : longblob   # Power at each dipole in the sourcemodel
ft    : longblob   # Full field trip results
%}


classdef Source< dj.Computed
    properties (Dependent)
        keySource
    end

    methods
        function v = get.keySource(~)
            % The epoch restriction is part of each SourceParm row, so it
            % cannot be applied once to the complete source relation.  Build
            % the union one parameter set at a time, retaining the SourceParm
            % key in each branch.
            base = proj(ns.Epoch * sloc.Leadfield * sloc.Headmodel * sloc.SourceParm);
            sourceParms = fetch(sloc.SourceParm, 'stag', 'parms');
            if isempty(sourceParms)
                v = base & 'FALSE';
                return
            end

            thisParm = sourceParms(1);
            epochRestriction = thisParm.parms.epoch;
            if ~isstruct(epochRestriction)
                epochRestriction = struct('etag', epochRestriction);
            end
            v = base & struct('stag', thisParm.stag) & epochRestriction;
            for iParm = 2:numel(sourceParms)
                thisParm = sourceParms(iParm);
                epochRestriction = thisParm.parms.epoch;
                if ~isstruct(epochRestriction)
                    epochRestriction = struct('etag', epochRestriction);
                end
                thisSource = base & struct('stag', thisParm.stag) & epochRestriction;
                v = base & (v | thisSource);
            end

        end
    end

    methods (Access = public)
        function GG = plotAverage(tbl,pv)
            arguments
                tbl
                pv.what (1,1) string = "3d"
                pv.equalClim (1,1) logical = true
                pv.fun =  @(x) mean(x,1,"omitmissing");
            end
            G= fetchtable(tbl,'power','label');
            GG= groupsummary(G,"label",pv.fun,"power");
            GG = sortrows(GG,"label");
            sourcemodel = fetch1(sloc.Headmodel & tbl,'sourcemodel','LIMIT 1');
            power =GG.fun1_power;
            ax = gobjects(height(GG),1);
            switch lower(pv.what)
                case "3d"
                    clf;
                    tiledlayout('flow')
                    for g= 1:height(GG)
                        ax(g) = nexttile;
                        scatter3(sourcemodel.pos(sourcemodel.inside,1), ...
                            sourcemodel.pos(sourcemodel.inside,2), ...
                            sourcemodel.pos(sourcemodel.inside,3), ...
                            35, ...
                            power(g,sourcemodel.inside),'filled');
                        axis equal
                        xlabel('x (mm)');
                        ylabel('y (mm)');
                        zlabel('z (mm)');
                        title (GG{g,"label"})
                        cb = colorbar;
                        ylabel(cb,'Source power');
                        view(3)
                        grid on
                        rotate3d on
                    end

                case "mri"
                    mri = fetch1(sloc.Headmodel & tbl,'mri','LIMIT 1');
                    sourceData = fetch1(tbl,'ft','LIMIT 1');
                    sourceData= sourceData.response;
                    for g = 1:height(GG)

                        cfg = [];
                        cfg.parameter = 'pow';
                        cfg.interpmethod = 'sphere_weighteddistance';
                        cfg.downsample = 2;
                        cfg.sphereradius = 20;
                        cfg.power = 2;
                        sourceData.pow = power(g,:)';
                        sourceData = rmfield(sourceData,setdiff(fieldnames(sourceData),["pos" "pow" "inside" "unit" "coordsys"]));
                        functional = ft_sourceinterpolate(cfg, sourceData, mri);

                        funColorLim = prctile(power(g,:),[1 99]);
                        cfg = [];
                        cfg.method = 'ortho';
                        cfg.funparameter = 'pow';
                        cfg.maskparameter = 'pow';
                        %cfg.opacitymap = funColorLim;
                        cfg.funcolorlim = funColorLim;
                        cfg.location = 'max';
                        cfg.interactive = 'yes';
                        figure('Name',sprintf('Average source power: %s', ...
                            GG{g,"label"}),'Color','w');
                        ft_sourceplot(cfg, functional);
                    end

                otherwise
                    error('Unknown plot type: %s',pv.what);
            end
            if ~isempty(ax)
                linkprop(ax,{"CameraPosition","CameraUpVector"});
                if pv.equalClim
                    cLim = cat(1,ax.CLim);
                    cLim = [min(cLim(:,1)) max(cLim(:,2))];
                    set(ax,'CLim',cLim)
                end
            end

        end
        function plot(tbl,pv)
            arguments
                tbl
                pv.what (1,:) string = "3d"
                pv.newFigureEach (1,:) = "subject"
                pv.equalClim (1,1) logical = false;
                pv.average (1,1) logical = false
            end
            current = dictionary;
            current(pv.newFigureEach) ="";
            ax = [];
            for tpl = fetch(tbl,'*',sprintf('ORDER BY %s',strjoin(pv.newFigureEach,',')))'
                for fld = pv.newFigureEach
                    if ~ismember(tpl.(fld),current(fld))
                        if ~isempty(ax)
                            if pv.equalClim
                                cLim = cat(1,ax.CLim);
                                cLim = [min(cLim(:,1)) max(cLim(:,2))];
                                set(ax,'CLim',cLim)
                            end

                            linkprop(ax,{"CameraPosition","CameraUpVector"});
                        end
                        hFig = figure;
                        ax  = [];
                        tiledlayout('flow');
                        current(fld) = tpl.(fld);
                        break;
                    end
                end
                sourcemodel =fetch1(sloc.Headmodel & tpl,'sourcemodel');

                if count(ns.DimensionCondition & tpl)>1;continue;end
                condition = fetch1(ns.DimensionCondition & tpl,'name');
                for wh = pv.what
                    switch lower(wh)
                        case "3d"
                            % Simple 3D dots scatter
                            ax  =[ax nexttile]; %#ok<AGROW>
                            scatter3(sourcemodel.pos(sourcemodel.inside,1), ...
                                sourcemodel.pos(sourcemodel.inside,2), ...
                                sourcemodel.pos(sourcemodel.inside,3), ...
                                35, ...
                                tpl.power(sourcemodel.inside),'filled');
                            axis equal
                            xlabel('x (mm)');
                            ylabel('y (mm)');
                            zlabel('z (mm)');
                            title (sprintf('%s epoch: %s (%s) in atlas %s',tpl.subject, tpl.etag,condition,tpl.hmtag));
                            cb = colorbar;
                            ylabel(cb,'Source power');
                            view(3)
                            grid on
                            rotate3d on
                        case "mri"
                            mri = fetch1(sloc.Headmodel & tpl,'mri');

                            cfg = [];
                            cfg.parameter = 'pow';
                            cfg.interpmethod ='sphere_weighteddistance';
                            cfg.downsample = 2;
                            cfg.sphereradius =20;
                            cfg.power = 2;
                            sourcedata = tpl.source;
                            sourcedata = rmfield(sourcedata,setdiff(fieldnames(sourcedata),["pos" "pow"  "unit" "coordsys"]));
                            functional = ft_sourceinterpolate(cfg,  ...
                                sourcedata,...
                                mri);

                            funColorLim = prctile(sourcedata.pow,[1 99]);

                            cfg = [];
                            cfg.method        = 'ortho';
                            cfg.funparameter  = 'pow';
                            % cfg.maskparameter = 'pow_abs';
                            cfg.funcolorlim   = funColorLim;
                            % cfg.opacitylim    = [0 colorScale];
                            % cfg.opacitymap    = 'rampup';
                            % cfg.location      = 'max';
                            cfg.interactive   = 'yes';
                            figure('Name','Source-localized power on MRI','Color','w');
                            ft_sourceplot(cfg, functional, mri);



                        otherwise
                            error("Unknown plot what :%s",wh);
                    end
                end
            end
            if ~isempty(ax)
                if pv.equalClim
                    cLim = cat(1,ax.CLim);
                    cLim = [min(cLim(:,1)) max(cLim(:,2))];
                    set(ax,'CLim',cLim)
                end
                linkprop(ax,{"CameraPosition","CameraUpVector"});
            end
        end
    end
    methods (Access = protected)
        function makeTuples(tbl,key)
            ft_defaults
            parms = fetch1(sloc.SourceParm & key, 'parms');
            % Get the data in Fieldtrip format
            data = ephys.fieldtrip.dataset(key);
            % Remove the reference if it is in the set
            cfg = [];
            cfg.channel = {'all','-E257'};
            data    = ft_selectdata(cfg, data);
            % Rereference to the average to match the FT headmodel
            cfg = [];
            cfg.reref = 'yes';
            cfg.refmethod = 'avg';
            cfg.refchannel = 'all';
            data = ft_preprocessing(cfg,data);

            % Pull the headmodel and leadfield
            [model] = fetch((sloc.Leadfield*sloc.Headmodel) & key,'leadfield','headmodel','sourcemodel','elec','mri','mesh');

            switch upper(parms.mode)
                case 'SNR'
                    % Use a single response window and determine power in
                    % the signal band for each source, divide it by signal
                    % in the noise band for each source. I.e., source
                    % localized SNR. 
                    % This uses DICS and 
                    for freq = fetch(ns.DimensionCondition & key,'*')'
                        f0 = freq.value;
                        if iscell(f0);f0=[f0{:}];end
                        f0 = unique(double(f0(:)'));
                        assert(isscalar(f0),"Only one frequency is allowed per condition.");

                        for harmonic = parms.harmonic(:)'
                            if isscalar(parms.harmonic)
                                label = string(f0);
                            else
                                label = string(f0) + "x" + string(harmonic);
                            end
                            f = f0*harmonic;

                            windowLength = diff(parms.responseWindow);
                            minimumF = 1000/windowLength;
                            if ~isfinite(f) || f < minimumF
                                fprintf(2, ...
                                    ['f0=%.6g Hz cannot be estimated in the shortest ' ...
                                    'analysis window (%.6g ms; minimum %.6g Hz).'], ...
                                    f, windowLength, minimumF);
                                continue;
                            end

                            frequencyBand = f + [-1 1]*max(parms.signalHalfWidth, parms.noiseHalfWidth);
                            signalBand= f + [-1 1]*parms.signalHalfWidth;
                            lowNoiseBand= f  -[parms.noiseHalfWidth parms.signalHalfWidth];
                            highNoiseBand= f   + [parms.signalHalfWidth parms.noiseHalfWidth];

                            assert(frequencyBand(1) >0, ...
                                "The frequency band must not extend below 0 Hz.");

                            cfg = [];
                            cfg.trials = find(ismember(data.trialinfo(:,1),freq.trials));
                            cfg.latency = parms.responseWindow/1000 + [0 -1./data.fsample];
                            dataResponse = ft_selectdata(cfg, data);

                            % Estimate all CSD bins needed for the signal and noise
                            % bands so that the common filter sees the same data.
                            cfg = parms.freq;
                            cfg.foilim = frequencyBand;
                            cfg.output = 'powandcsd';
                            freqResponse = ft_freqanalysis(cfg, dataResponse);

                            signalMask = do.ifinrange(freqResponse.freq,signalBand);
                            noiseHighMask = do.ifinrange(freqResponse.freq ,highNoiseBand);
                            noiseLowMask = do.ifinrange(freqResponse.freq ,lowNoiseBand);
                            assert(any(signalMask) && any(noiseLowMask) && any(noiseHighMask), ...
                                "The requested signal and noise widths do not contain frequency bins.");

                            cfg = parms.source;
                            cfg.frequency = frequencyBand;
                            cfg.sourcemodel = model.sourcemodel;
                            cfg.sourcemodel.leadfield = model.leadfield.leadfield;
                            cfg.sourcemodel.label = model.leadfield.label;
                            cfg.headmodel = model.headmodel;
                            source = ft_sourceanalysis(cfg, freqResponse);

                            % Set the common filters
                            cfg.sourcemodel.filter = source.avg.filter;
                            cfg.sourcemodel.filterdimord = source.avg.filterdimord;
                            cfg.sourcemodel.label = source.avg.label;
                            % And apply to each of the frequency ranges

                            cfg.frequency       = frequencyRange(freqResponse.freq(signalMask));
                            sourceSignal        = ft_sourceanalysis(cfg, freqResponse);
                            cfg.frequency = frequencyRange(freqResponse.freq(noiseLowMask));
                            sourceNoiseLow = ft_sourceanalysis(cfg, freqResponse);
                            cfg.frequency = frequencyRange(freqResponse.freq(noiseHighMask));
                            sourceNoiseHigh = ft_sourceanalysis(cfg, freqResponse);

                            nNoiseLow = nnz(noiseLowMask);
                            nNoiseHigh = nnz(noiseHighMask);
                            sourceNoise = (nNoiseLow*sourceNoiseLow.avg.pow + ...
                                nNoiseHigh*sourceNoiseHigh.avg.pow)/(nNoiseLow+nNoiseHigh);

                            snr = sourceSignal.avg.pow./sourceNoise;

                            tpl = key;
                            tpl.power = snr';
                            tpl.label = label;
                            tpl.ft = struct('response',sourceSignal, ...
                                'noiseLow',sourceNoiseLow,'noiseHigh',sourceNoiseHigh, ...
                                'common',source);
                            insert(tbl,makeMymSafe(tpl));
                        end
                    end

                case 'BASELINE'
                    % DICS analysis using SSVEPs wiht a common filter
                    % applied to a window before stimulus onset (baseWindow) and
                    % a window after stimulus onset (responseWindow).
                    % The power is the ratio of the stimulus power to
                    % base window power. 
                    % 
                    % Compared to the SNR mode, this mode normalizes to the
                    % same (signal) frequency in a different time window,
                    % not to the different (noise) frequencies in the same
                    % time window. 
                    %
                    % The power in the table is the contrast between the two windows:
                    % (respons-base)/base.
                    %
                    % The epochs should be determined for different SSVEP frequencies; the
                    % value of the ns.DimensionCondition associated with the epoch should be the SSVEP frequency.
                    %
                    % Frequencies that are too low to be  estimated in the base/response windows will be skipped.
                    %
                    % The dimension should correspond to the different
                    % frequencies; we loop over those condition (frequencies) here
                    %
                    for freq = fetch(ns.DimensionCondition & key,'*')'

                        f0 = freq.value;
                        if iscell(f0);f0=[f0{:}];end
                        f0 = unique(double(f0(:)'));
                        assert(isscalar(f0),"Only one frequency is allowed per condition.");
                        if ~iscell(parms.responseWindow)
                            parms.responseWindow = {parms.responseWindow};
                        end
                        assert(all(diff(parms.baseWindow)==diff(cat(1,parms.responseWindow{:}),1,2)), ...
                            "Base and response windows must have the same length to use a common filter.");


                        for harmonic = parms.harmonic(:)'
                            if isscalar(parms.harmonic)
                                label = string(f0);
                            else
                                label = string(f0) + "x" + string(harmonic);
                            end
                            f = f0*harmonic;

                            windowLength = diff(parms.baseWindow);
                            minimumF = 1000/windowLength;
                            if ~isfinite(f) || f < minimumF
                                fprintf(2, ...
                                    ['f0=%.6g Hz cannot be estimated in the shortest ' ...
                                    'analysis window (%.6g ms; minimum %.6g Hz).'], ...
                                    f, windowLength, minimumF);
                                continue;
                            end

                            signalBand= f + [-1 1]*parms.signalHalfWidth;
                            assert(signalBand(1) >0, ...
                                "The frequency band must not extend below 0 Hz.");

                            cfg = [];
                            cfg.trials = find(ismember(data.trialinfo(:,1),freq.trials));
                            cfg.latency = parms.baseWindow/1000 + [0 -1./data.fsample];
                            dataBase = ft_selectdata(cfg, data);

                            dataResponse = cell(numel(parms.responseWindow),1);
                            for r = 1:numel(parms.responseWindow)
                                cfg.latency = parms.responseWindow{r}/1000 + [0 -1./data.fsample];
                                dataResponse{r} = ft_selectdata(cfg, data);
                            end

                            % Combine the base window with one response window
                            % to estimate the filter.
                            cfg = [];
                            dataCommon = ft_appenddata(cfg, dataBase, dataResponse{1});
                            % And combine all response windows to apply the
                            % filter (And determine power)
                            dataResponse = ft_appenddata(cfg, dataResponse{:});

                            % Estimate all CSD bins needed for the signal and noise
                            % bands so that the common filter sees the same data.
                            cfg = parms.freq;
                            cfg.foilim = signalBand;
                            cfg.output = 'powandcsd';
                            freqCommon = ft_freqanalysis(cfg, dataCommon);
                            freqResponse = ft_freqanalysis(cfg, dataResponse);
                            freqBase  = ft_freqanalysis(cfg, dataBase);

                            signalMask = do.ifinrange(freqResponse.freq,signalBand);
                            assert(any(signalMask), "The requested signal widths do not contain frequency bins.");

                            cfg = parms.source;
                            cfg.frequency = signalBand;
                            cfg.sourcemodel = model.sourcemodel;
                            cfg.sourcemodel.leadfield = model.leadfield.leadfield;
                            cfg.sourcemodel.label = model.leadfield.label;
                            cfg.headmodel = model.headmodel;
                            source = ft_sourceanalysis(cfg, freqCommon);

                            % Set the common filters
                            cfg.sourcemodel.filter = source.avg.filter;
                            cfg.sourcemodel.filterdimord = source.avg.filterdimord;
                            cfg.sourcemodel.label = source.avg.label;
                            % And apply to each of the frequency ranges in the
                            % signal and base windows

                            cfg.frequency       = frequencyRange(freqResponse.freq(signalMask));
                            sourceSignal        = ft_sourceanalysis(cfg, freqResponse);
                            sourceSignalBase    = ft_sourceanalysis(cfg, freqBase);
                        
                            power= sourceSignal.avg.pow./sourceSignalBase.avg.pow;

                            tpl = key;
                            tpl.power = power';
                            tpl.label = label;
                            tpl.ft = struct('response',sourceSignal, ...                                
                                'common',source,...
                                'base',sourceSignalBase);
                            insert(tbl,makeMymSafe(tpl));
                        end
                    end

                otherwise
                    error('Unknown source localization mode %s',parms.mode);
            end

        end
    end


end

function frequency = frequencyRange(frequencies)
% Return the compact FieldTrip frequency selector for a set of bins.
frequencies = frequencies(:)';
if isscalar(frequencies)
    frequency = frequencies;
else
    frequency = [frequencies(1) frequencies(end)];
end
end
