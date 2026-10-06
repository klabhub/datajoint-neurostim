%{
# Leadfields for source localization in one session.
-> sloc.Headmodel  
-> ns.Session    
-> sloc.LeadfieldParm
---
leadfield : longblob # The leadfield in FieldTrip format
elec : longblob # The electrode positions used in the leadfield computation, in FieldTrip format
movement :float # The median movement of the electrodes after projecting to the scalp (in mm)
%}
% The headmodel is computed once per subject and per electrode montage (i.e. for a generic 256 channel net)
% The leadfield depends on the locations of the electrodes, hence they differ per session.
% The elec are projected to the scalp, and then snapped to the nearest surface vertex for the leadfield computation.

classdef Leadfield < dj.Computed & dj.DJInstance
    properties (Dependent)
        keySource
    end

    methods
        function v= get.keySource(~)
            % Restrict to sessions with an experiment in one of the paradigms listed for each parameter set.
            % proj() reduces to the primary key so each combination appears once.
            % The restriction must be applied after proj; a proj on top of a restricted
            % relation generates invalid SQL (two WHERE clauses) when populate adds its own.
            paradigms = proj(sloc.LeadfieldParmParadigm,'name->paradigm');
            v = proj(sloc.Headmodel * ns.Session * sloc.LeadfieldParm) & (ns.Experiment * paradigms);
        end

        function plot(tbl,pv)
            arguments
                tbl (1,1) sloc.Leadfield
                pv.what (1,:) string ="norm"
            end

            [leadfield, elec] = fetch1(tbl, 'leadfield', 'elec');
            mesh = fetch1(sloc.Headmodel & tbl, 'mesh');

            for what = pv.what
                switch what
                    case "norm"
                        % Show leadfield locations and strength in relation
                        % to the anatomical mesh and recording electrodes.
                        ft_plot_mesh(mesh, 'edgecolor', 'none', ...
                            'facecolor', 'skin_medium_light', ...
                            'facealpha', 0.25, 'surfaceonly', true);
                        hold on

                        inside = sourceInside(leadfield);
                        strength = nan(size(leadfield.pos, 1), 1);
                        for dipole = find(inside)'
                            if ~isempty(leadfield.leadfield{dipole})
                                strength(dipole) = norm(leadfield.leadfield{dipole});
                            end
                        end
                        valid = inside & isfinite(strength) & strength > 0;
                        color = strength(valid);
                        scatter3(leadfield.pos(valid,1), leadfield.pos(valid,2), ...
                            leadfield.pos(valid,3), 24, color, 'filled');

                        clim([min(color) max(color)])
                        colormap hot
                        nElectrodes = 0;
                        if isfield(elec, 'elecpos')
                            nElectrodes = size(elec.elecpos, 1);
                            plot3(elec.elecpos(:,1), elec.elecpos(:,2), ...
                                elec.elecpos(:,3), 'k.', 'MarkerSize', 10);
                        end
                        ft_plot_axes(mesh)
                        axis equal tight vis3d
                        view(3)
                        cb = colorbar;
                        cb.Label.String = 'leadfield norm';
                        title(sprintf('Leadfield strength (%d dipoles, %d electrodes)', ...
                            nnz(valid), nElectrodes))
                        hold off
                    otherwise
                        error('Unknown plot type: %s.\n', what);
                end
            end
        end
    end 
    methods (Access=protected)


        function makeTuples(tbl,key)
            assert(exist("ft_defaults","file"),"Please add FieldTrip to your path.")
            if ~exist("ft_version","file"  )
                ft_defaults;                            
            end
            ft_hastoolbox('simbio', 1);   % puts external/simbio (sb_find_elec, sb_calc_vecx, sb_solve) on the path
          
            fidLabels = {'nas','lpa','rpa'}; % CTF only for now
            noLeadfieldLabels = cat(2,fidLabels,{'Cz'}); % Cz is reference, no leadfield needed
           
            [mri, headmodel, scalp,sourcemodel] = fetch1(sloc.Headmodel & key, 'mri', 'headmodel', 'scalp','sourcemodel');
            assert(isfield(mri,'cfg') && isfield(mri.cfg,'fiducial'),'This MRI does not specify the fiducials.');
            parms  = fetch1(sloc.LeadfieldParm & key, 'parms');

            gpsFile = ns.File & 'extension=".gpsr"' & 'filename LIKE "%solved%"' & key;
            if ~exists(gpsFile)
                % Maybe it exists but has not been added  to the ns.File
                % table.
                f = fullfile(folder(ns.Session & key),key.subject + "*.solved.gpsr");
                solvedGps = dir(f);
                if isempty(solvedGps)
                    % Test to see whetehr there is a GPS file at all
                    f = fullfile(folder(ns.Session & key),key.subject + "*.gpsr");
                    unsolvedGps = dir(f);
                    if isempty(unsolvedGps)
                        error('No solved GPS file found for this session. (%s/%s)',key.subject,key.session_date);
                    else
                        error("Session %s for %s has an unsolved GPS file. Solve it, then retry the leadfield computation\n",key.subject,key.session_date);
                    end
                else
                    % Add it, then requery
                    updateWithFiles(ns.File,key,solvedGps);                    
                    gpsFile = ns.File & 'extension=".gpsr"' & 'filename LIKE "%solved%"' & key;            
                end
            end

            
            
            filename= fullfile(folder(ns.Session & key),gpsFile{1,"filename"},'coordinates.gpsc');
            assert(exist(filename,"file"), '%s not found',filename);
            GPS= readtable(filename, FileType="text");
            GPS = renamevars(GPS,"Var"  + string(1:4),["label" "Y" "X" "Z"]);
            % FT wants Y X Z
            elec.chanpos = [GPS.Y GPS.X GPS.Z];
            elec.chantype = repmat({'eeg'},height(GPS),1);
            elec.chanunit = repmat({'cm'},height(GPS),1);
            elec.elecpos = elec.chanpos;
            elec.pnt = elec.chanpos;
            elec.label   = GPS.label;
             % Rename EGI fiducials to FieldTrip-compatible names
            elec.label(strcmpi(elec.label, 'FidNz'))  = {'nas'};
            elec.label(strcmpi(elec.label, 'FidT9'))  = {'lpa'};
            elec.label(strcmpi(elec.label, 'FidT10')) = {'rpa'};
            elec = ft_convert_units(elec,'mm');
            elec = ft_datatype_sens(elec);
            
      
            %% Align the fiducials in the net with the fiducials in the MRI.
            cfg = [];
            elec.fid.label = fidLabels;
            cfg.target.label =fidLabels;

            fid_vox  = [];   
            for lbl = string(fidLabels)
                % Find the fiducial label in the elec then assign the pos
                elec.fid.pos = elec.elecpos(strcmpi(lbl,elec.label),:);
                % Get the fiducial position from the MRI
                fid_vox  = [fid_vox; mri.cfg.fiducial.(lbl)]; %#ok<AGROW>
            end
            % Transform the mri coords to ctf
            fid_ctf = ft_warp_apply(mri.transform, fid_vox, 'homogeneous');
            cfg.target.chanpos = fid_ctf;
            cfg.target.elecpos = fid_ctf;
            cfg.target.pnt      = fid_ctf;
            cfg.target.unit     = 'mm';
            cfg.method = 'fiducial';
            cfg.elec  =elec;
            elec = ft_electroderealign(cfg);

            %% Fine tune by projecting the electrodes to the scalp
            cfg = [];
            cfg.method    = 'project';
            cfg.headshape = scalp;
            elec_projected = ft_electroderealign(cfg,elec);
            % Sanity check
            movement = sqrt(sum((elec_projected.elecpos - elec.elecpos).^2, 2));
            fprintf('Scalp projection median movement: %.2f mm (Quintiles: %.2f, %.2f, %.2f, %.2f, %.2f)\n', median(movement), quantile(movement, 0.2), quantile(movement, 0.4), quantile(movement, 0.6), quantile(movement, 0.8), quantile(movement, 1.0    ));            
            
            %% Prune to keep only the real electrodes
             % Keep only the real electrodes (no fiducials, no reference)
            keep = ~ismember(elec_projected.label, noLeadfieldLabels);
            for f = ["label" "elecpos" "chanpos" "chantype" "chanunit" "pnt"]
                if isfield(elec_projected,f)
                    elec_projected.(f) = elec_projected.(f)(keep,:);
                end
            end
            % Snap the electrodes to the head model
            elec_projected = snapElectrodes(elec_projected, headmodel);
            nrElectrodes  = numel(elec_projected.label);
            %% Compute the leadfield
            pool = nsParPool;                   % [] => serial , otherwise parallel on the current pool
            
            if strcmpi(parms.mode, "auto")
                % Determine whether e2d or d2e is more efficient.  Only source
                % positions marked inside are actually solved in d2e mode.
                if isfield(sourcemodel, 'inside')
                    inside = sourcemodel.inside;
                    if islogical(inside)
                        nDipoles = nnz(inside);
                    else
                        nDipoles = numel(inside);
                    end
                else
                    nDipoles = size(sourcemodel.pos, 1);
                end
                if nrElectrodes > nDipoles
                    parms.mode = 'd2e'; % More efficient to solve from dipole positions
                else
                    parms.mode = 'e2d'; % More efficient to solve from electrode positions
                end
            end
            % Call the appropriate leadfield computation function
            % Note that each function returns identical results, they differ only in computational efficiency.
            switch parms.mode
                case 'ft'
                    % Use the FieldTrip method (no parallelization, always per electrode)                
                    cfg = parms.cfg;
                    cfg.channel     = elec_projected.label;
                    cfg.elec        = elec_projected;
                    cfg.headmodel   = headmodel;
                    cfg.sourcemodel = sourcemodel;
                    leadfield = ft_prepare_leadfield(cfg);  
                case 'e2d'
                    % Electrode-to-dipole
                    % One FE solve per electrode; more efficient when the sourcemodel is large.
                    cfg = parms.cfg;
                    cfg.channel     = elec_projected.label;
                    cfg.elec        = elec_projected;
                    cfg.headmodel   = headmodel;
                    cfg.sourcemodel = sourcemodel;
                    leadfield = leadfieldPerElectrode(cfg, pool);
                case 'd2e'
                    % Dipole-to-electrode     
                    % One FE solve per dipole orientation (3 per source position); much cheaper than one per
                    % electrode when the sourcemodel is small.                    
                    cfg = parms.cfg;
                    cfg.channel     = elec_projected.label;
                    cfg.elec        = elec_projected;
                    cfg.headmodel   = headmodel;
                    cfg.sourcemodel = sourcemodel;
                    leadfield = leadfieldPerDipole(cfg, pool);
            end            

            tpl = key;
            tpl.movement = median(movement);
            tpl.leadfield = leadfield;
            tpl.elec = elec_projected;
            insert(tbl, tpl);

        end
    end
end

function inside = sourceInside(sourcemodel)
if isfield(sourcemodel, 'inside')
    inside = sourcemodel.inside(:);
    if ~islogical(inside)
        mask = false(size(sourcemodel.pos, 1), 1);
        mask(inside) = true;
        inside = mask;
    end
else
    inside = true(size(sourcemodel.pos, 1), 1);
end
end

function elec = snapElectrodes(elec, headmodel)
% SimBio needs the electrodes on mesh nodes: snap to the nearest vertex of the outer surface
% (same snapping as ft_prepare_vol_sens)
surfNodes = surfaceNodes(headmodel);
surfPos   = headmodel.pos(surfNodes,:);
for j = 1:size(elec.elecpos,1)
    [~, k] = min(sum((surfPos - elec.elecpos(j,:)).^2, 2));
    elec.elecpos(j,:) = surfPos(k,:);
end
elec.chanpos = elec.elecpos;
if isfield(elec,'pnt'), elec.pnt = elec.elecpos; end
end

function nodes = surfaceNodes(hm)
% Vertices of the faces that belong to exactly one element (as in FieldTrip's private mesh2edge).
if isfield(hm,'hex')
    e = hm.hex;
    f = cat(1, e(:,[1 2 3 4]), e(:,[5 6 7 8]), e(:,[1 2 6 5]), ...
        e(:,[2 3 7 6]), e(:,[3 4 8 7]), e(:,[4 1 5 8]));
elseif isfield(hm,'tet')
    e = hm.tet;
    f = cat(1, e(:,[1 2 3]), e(:,[2 3 4]), e(:,[3 4 1]), e(:,[4 1 2]));
else
    error('Unsupported SimBio mesh: expected hex or tet elements.');
end
[~, ~, ic] = unique(sort(f,2), 'rows');
single = accumarray(ic, 1) == 1;
nodes  = unique(reshape(f(single(ic),:), [], 1));
end

function row = solveElectrode(stiff, elecnode, refNode)
% Potential on all mesh nodes for a unit source at elecnode, relative to refNode.
% Runs on the client or on a pool worker.
if isa(stiff,'parallel.pool.Constant'), stiff = stiff.Value; end

vecb = zeros(size(stiff,1),1);
vecb(elecnode) = 1;
row = reshape(sb_calc_vecx(stiff, vecb, refNode), 1, []);
end

function reportProgress(n, total, t0)
el  = duration(seconds(toc(t0)));
eta = duration(el/n*(total-n),'Format','hh:mm:ss');
fprintf('FE solves: %d/%d (%.0f%%), elapsed %s , ETA %s \n', n, total, 100*n/total, el, eta);
end


%%
% Functions for computing the leadfield per dipole.

function sourcemodel = leadfieldPerDipole(cfg, pool)
% Leadfield for a SimBio FE headmodel by solving one FE system per dipole orientation.
%
% This is an alternative to ft_prepare_leadfield with a precomputed transfer matrix (one FE
% solve per electrode). With few source positions (3 solves per dipole) it needs far fewer solves than
% electrodes. The stiffness matrix is scaled and preconditioned once, and then reused for all solves.
% The FE system is symmetric, so the potential at the electrodes for a dipole is identical to
% transfer(elec,:)*rhs, i.e. the result equals that of ft_prepare_leadfield (up to solver tolerance).
%
% sourcemodel = sloc.leadfieldPerDipole(cfg, pool)
%
% cfg fields (same as the cfg for ft_prepare_leadfield):
%   cfg.headmodel    SimBio headmodel with pos, hex or tet, and stiff
%   cfg.elec         electrodes (already on the scalp), the first is the reference
%   cfg.channel      channel labels (subset of elec.label), in the order of the output
%   cfg.sourcemodel  pos (N*3) and optionally inside
%   cfg.normalize, cfg.normalizeparam, cfg.weight, cfg.reducerank, cfg.backproject
% pool (optional) parallel pool; solves run on it if it has workers, otherwise serially.

if nargin < 2, pool = []; end
ft_hastoolbox('simbio', 1);

headmodel   = cfg.headmodel;
sourcemodel = cfg.sourcemodel;
elec        = cfg.elec;

[found, chanIdx] = ismember(cfg.channel, elec.label);
assert(all(found), 'Not all channels are in cfg.elec.');
assert(chanIdx(1) == 1, 'The first channel must be the reference electrode (first label in cfg.elec).');

elecnodes = sb_find_elec(headmodel, elec);
elecnodes = elecnodes(:);
refNode   = elecnodes(1);
elecnodes = elecnodes(chanIdx);

if isfield(sourcemodel, 'inside')
    inside = sourcemodel.inside;
    if ~islogical(inside)
        tmp = false(size(sourcemodel.pos,1),1);
        tmp(inside) = true;
        inside = tmp;
    end
else
    inside = true(size(sourcemodel.pos,1),1);
end
inside   = inside(:);
insideIx = find(inside);
nDip     = numel(insideIx);

opt.reducerank     = getopt(cfg, 'reducerank', 'no');
opt.backproject    = getopt(cfg, 'backproject', 'yes');
opt.normalize      = getopt(cfg, 'normalize', 'no');
opt.normalizeparam = getopt(cfg, 'normalizeparam', 0.5);
weight             = getopt(cfg, 'weight', []);
if isscalar(weight), weight = weight*ones(nDip,1); end

fprintf('Per-dipole leadfield: %d dipoles, %d channels, %d mesh nodes\n', nDip, numel(elecnodes), size(headmodel.pos,1));
tStart = tic;
solver = prepareSolver(headmodel.stiff, refNode);
fprintf('Solver prepared in %.0f s\n', toc(tStart));

% Venant load vectors (n nodes x 3 orientations) for each dipole
dirs = eye(3);
rhs  = cell(nDip,1);
for i = 1:nDip
    rhs{i} = sb_rhs_venant(repmat(sourcemodel.pos(insideIx(i),:),3,1), dirs, headmodel);
end

pot = cell(nDip,1);     % each nElec x 3
if ~isempty(pool) 
    solverConst = parallel.pool.Constant(solver);
    futures(nDip,1) = parallel.FevalFuture;
    for i = 1:nDip
        futures(i) = parfeval(pool, @solveDipole, 1, solverConst, rhs{i}, elecnodes);
    end
    cancelFutures = onCleanup(@() cancel(futures)); 
    for n = 1:nDip
        [i, p] = fetchNext(futures);
        pot{i} = p;
        reportProgress(n, nDip, tStart);
    end
else
    for i = 1:nDip
        pot{i} = solveDipole(solver, rhs{i}, elecnodes);
        reportProgress(i, nDip, tStart);
    end
end

sourcemodel.leadfield = cell(1, size(sourcemodel.pos,1));
sourcemodel.leadfield(:) = {[]};
for i = 1:nDip
    w = 1;
    if ~isempty(weight), w = weight(i); end
    sourcemodel.leadfield{insideIx(i)} = postprocess(pot{i}, opt, w);
end
sourcemodel.inside          = inside;
sourcemodel.label           = cfg.channel;
sourcemodel.leadfielddimord = '{pos}_chan_ori';
end

function sourcemodel = leadfieldPerElectrode(cfg, pool)
% Leadfield for a SimBio FE headmodel by solving one FE system per electrode.
%
% This is the transfer-matrix route used by ft_prepare_leadfield. The FE
% system is solved once for each non-reference electrode, after which
% ft_prepare_leadfield applies the transfer matrix to the source positiosloc.
%
% sourcemodel = sloc.leadfieldPerElectrode(cfg, pool)
%
% cfg fields are the same as for ft_prepare_leadfield. pool (optional) is a
% parallel pool; solves run on it if it has workers, otherwise serially.

if nargin < 2, pool = []; end
ft_hastoolbox('simbio', 1);

headmodel = cfg.headmodel;
elec      = cfg.elec;
nrLabels  = numel(elec.label);

elecnodes = sb_find_elec(headmodel, elec);

% One FE solve per electrode label (reference = first label, as in sb_transfer).
stiff       = headmodel.stiff;
nPos        = size(headmodel.pos, 1);
T           = zeros(nrLabels, nPos);
solveLabels = 2:nrLabels;
nSolve      = numel(solveLabels);
refNode     = elecnodes(1);
tStart      = tic;

% Progress is printed from the client loop because DataQueue callback output
% does not reach the terminal.
if isempty(pool)
    for k = 1:nSolve
        T(solveLabels(k),:) = solveElectrode(stiff, elecnodes(solveLabels(k)), refNode);
        reportProgress(k, nSolve, tStart);
    end
else
    stiffConst = parallel.pool.Constant(stiff);
    futures(nSolve,1) = parallel.FevalFuture;
    for k = 1:nSolve
        futures(k) = parfeval(pool, @solveElectrode, 1, stiffConst, elecnodes(solveLabels(k)), refNode);
    end
    cancelFutures = onCleanup(@() cancel(futures)); 
    for n = 1:nSolve
        [k, row] = fetchNext(futures);
        T(solveLabels(k),:) = row;
        reportProgress(n, nSolve, tStart);
    end
end

headmodel.transfer = T;
headmodel.elec     = elec;
cfg.headmodel      = headmodel;
sourcemodel        = ft_prepare_leadfield(cfg);
end

function lf = postprocess(lf, opt, w)
% Same steps as ft_compute_leadfield for a single dipole (EEG, average reference).
lf = lf - mean(lf,1);

switch opt.reducerank
    case 'yes', r = 2;
    case 'no',  r = 3;
    otherwise,  r = opt.reducerank;
end
if r < 3
    [u, s, v] = svd(lf);
    d = diag(s);
    s(:) = 0;
    for j = 1:r, s(j,j) = d(j); end
    if istrue(opt.backproject)
        lf = u*s*v';
    else
        lf = lf*v(:,1:r);
    end
end

switch opt.normalize
    case 'yes'
        if opt.normalizeparam == 0.5
            nrm = norm(lf,'fro');
        else
            nrm = sum(lf(:).^2)^opt.normalizeparam;
        end
        if nrm > 0, lf = lf./nrm; end
    case 'column'
        for j = 1:size(lf,2)
            lf(:,j) = lf(:,j)./(sum(lf(:,j).^2)^opt.normalizeparam);
        end
end
lf = lf*w;
end

function v = getopt(cfg, name, default)
if isfield(cfg, name) && ~isempty(cfg.(name))
    v = cfg.(name);
else
    v = default;
end
end

function S = prepareSolver(stiff, refNode)
% One-time setup for the solves, following sb_calc_vecx and sb_solve: Dirichlet condition (zero potential)
% at the reference node, diagonal scaling and incomplete Cholesky preconditioner.
n = size(stiff,1);
vecdi = zeros(n,1);
vecdi(refNode) = 1;
stiff = sb_set_bndcon(stiff, zeros(n,1), vecdi, zeros(n,1));
dkond = 1./sqrt(diag(stiff));
[ii, jj, s] = find(stiff);
clear stiff
s = (s.*dkond(ii)).*dkond(jj);
s(1) = 1;
L = sparse(ii, jj, s, n, n, length(s));
try
    L = ichol(L);
catch
    disp('Could not compute incomplete Cholesky-decomposition. Rescaling stiffness matrix...')
    alpha = 1/(0.5e-6*8 + 1);
    s = alpha*s;
    dia = ii == jj;
    s(dia) = sqrt(s(dia)/alpha);
    s(1) = 1;
    L = ichol(sparse(ii, jj, s, n, n, length(s)));
end
A = sparse(ii, jj, s, n, n, length(s));
clear ii jj s
A = A + A' - sparse(1:n, 1:n, diag(A), n, n, n);
S.A = A;
S.L = L;
S.Lt = L';
S.dkond = dkond;
S.refNode = refNode;
end

function pot = solveDipole(S, rhs, elecnodes)
% Potentials at the electrode nodes for the 3 orientations of one dipole.
if isa(S,'parallel.pool.Constant'), S = S.Value; end
pot = zeros(numel(elecnodes), size(rhs,2));
for k = 1:size(rhs,2)
    b = full(rhs(:,k));
    b(S.refNode) = 0;
    b = b.*S.dkond;
    x0 = S.Lt \ (S.L \ (-b));
    [x, flag, relres] = pcg(S.A, b, 10e-9, 5000, S.L, S.Lt, x0);
    if flag ~= 0
        warning('leadfieldPerDipole:pcg', 'pcg did not converge (flag %d, relres %.2g).', flag, relres);
    end
    x = x.*S.dkond;
    pot(:,k) = x(elecnodes);
end
end

