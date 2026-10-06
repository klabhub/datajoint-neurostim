%{
# Headmodel for source localization
->sloc.Fiducial
->sloc.HeadmodelParm
---
mri : longblob # The MRI volume in FieldTrip format
segmented: longblob # The segmented MRI in FieldTrip format
mesh : longblob # The mesh in FieldTrip format
headmodel : longblob # The headmodel in FieldTrip format
scalp : longblob # The scalp mesh in FieldTrip format
sourcemodel : longblob # The sourcemodel in FieldTrip format
%}

classdef Headmodel < dj.Computed & dj.DJInstance
    properties (Dependent   )
        atlas
    end

    methods
        function tissueLabels = get.atlas(tbl)
            for key = fetch(tbl * sloc.HeadmodelParm,'parms')'
                if key.parms.sourcemodel.cfg.method == "basedonmni"
                    [~,ftpath] = ft_version;
                    atlasFile = fullfile(ftpath, 'template/atlas/', key.parms.sourcemodel.atlas);                    
                    assert(exist(atlasFile,"file"), 'Atlas file (%s) not found', atlasFile);
                    thisAtlas = ft_read_atlas(atlasFile);
                    tissueLabels = thisAtlas.tissuelabel;
                end
            end
        end
    end
    methods(Access=public)
        function plot(tbl,pv)
            arguments
                tbl (1,1) sloc.Headmodel
                pv.what (1,:) string 
            end

            [mri, mesh, sourcemodel] = fetch1(tbl, 'mri', 'mesh', 'sourcemodel');
            parms = fetch1(sloc.HeadmodelParm & tbl, 'parms');
            
            for what = pv.what
                switch what
                    case "segmentation"
                        % convert from probabilistic/binary into indexed representation
                        segmentedmri_indexed = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
                        % also add the anatomical mri
                        segmentedmri_indexed.anatomy = mri.anatomy;
                        cfg              = [];
                        cfg.anaparameter = 'anatomy';
                        cfg.funparameter = 'tissue';
                        cfg.funcolormap  = lines(6);              % distinct color per tissue + background
                        cfg.location     = 'center';
                        ft_sourceplot(cfg, segmentedmri_indexed);
                    case "mesh"
                        ft_plot_mesh(mesh, 'edgecolor','none', 'facecolor', 'skin_medium_light', 'facealpha', 0.7,'surfaceonly',true);
                        ft_plot_axes(mesh)
                        alpha 1
                        material default
                        camlight
                    case "allmeshes"
                        tiledlayout('flow')
                        for tissue =string(mesh.tissuelabel)'
                            nexttile
                            tissueMesh = rmfield(mesh,{'tissue','tissuelabel'});
                            tissueMesh.hex = mesh.hex(mesh.tissue==find(strcmp(tissue,mesh.tissuelabel)),:);
                            ft_plot_mesh(tissueMesh, 'edgecolor','none', 'facealpha', 0.7,'surfaceonly',false);
                            ft_plot_axes(tissueMesh)
                            title (tissue)
                            alpha 1
                            material default
                            camlight
                        end
                    case "sourcemodel"
                        % Show the dipoles retained by the source model in
                        % relation to the anatomical mesh.
                        inside = sourcemodel.inside(:);
                        pos = sourcemodel.pos(inside, :);
                        hold on
                        if isfield(mri, 'anatomy') && isfield(mri, 'transform')
                            sliceLocation = mean(pos, 1);
                            anat = double(mri.anatomy);
                            % Robust contrast: clip the bright outliers that otherwise
                            % make the whole volume look uniformly dark.
                            clim = prctile(anat(:), [1 99]);
                            for orientation = eye(3)
                                ft_plot_slice(anat, ...
                                    'transform', mri.transform, ...
                                    'location', sliceLocation, ...
                                    'orientation', orientation', ...
                                    'colormap', gray(256), ...
                                    'clim', clim, ...
                                    'interpmethod', 'linear', ...
                                    'resolution', 1, ...
                                    'doscale', false, ...
                                    'facealpha', 1,...
                                    'unit','mm');
                            end
                        end
                        ft_plot_mesh(mesh, 'edgecolor', 'none', ...
                            'facecolor', 'skin_medium_light', ...
                            'facealpha', 0.25, 'surfaceonly', true);
                        plot3(pos(:,1), pos(:,2), pos(:,3), '.', ...
                            'Color', [0.8500 0.3250 0.0980], ...
                            'MarkerSize', 10);
                        if string(parms.sourcemodel.cfg.method) == "basedonmni"
                            tissueLabels = string(tbl.atlas);
                            if isfield(sourcemodel, 'tissuelabel')
                                tissueLabels = string(sourcemodel.tissuelabel);
                            end
                            assert(numel(tissueLabels) == size(pos, 1), ...
                                'MNI sourcemodel has %d dipoles but %d atlas tissue labels.', ...
                                size(pos, 1), numel(tissueLabels));
                            text(pos(:,1), pos(:,2), pos(:,3), tissueLabels(:), ...
                                'FontSize', 8, 'Color', [0 0 0], ...
                                'Interpreter', 'none', 'VerticalAlignment', 'bottom');
                        end
                        ft_plot_axes(mesh)
                        axis equal tight vis3d
                        view(3)
                        title(sprintf('Sourcemodel dipoles inside mesh (%d)', size(pos, 1)))
                        hold off
                        legend off
                    otherwise
                    error( 'Unknown plot type: %s.\n', what);
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
            if ispc
                fprintf(2,"The fortran module that creates the simbio headmodel is known to be buggy on Windows. Probably best to run this on linux.\n ");
            else
                fprintf(2,"Assuming Fortran libraries are available on this system.\n ");
                % On our HPC system something like this is needed before starting matlab:
                % module load gcc/4.9.4 && GFORTRAN_LIB=$(dirname "$(gcc -print-file-name=libgfortran.so.3)") && export LD_LIBRARY_PATH="$GFORTRAN_LIB:$LD_LIBRARY_PATH"
            end
           
            parms = fetch1(sloc.HeadmodelParm & key,'parms');
            fids = fetch(sloc.Fiducial & key,'*');


            %% Check that we can actually run till the end.
            [~, ftpath] = ft_version;
            if (lower(parms.sourcemodel.cfg.method)=="basedonmni")
                % Needs template and atlas files
                templateFile = fullfile(ftpath, 'template/sourcemodel/', parms.sourcemodel.template);
                if ~endsWith(templateFile,'.mat');templateFile = [templateFile  '.mat'];end
                atlasFile = fullfile(ftpath, 'template/atlas/', parms.sourcemodel.atlas);
                assert(exist(templateFile,"file"), 'Template file (%s) not found', templateFile);
                assert(exist(atlasFile,"file"), 'Atlas file (%s) not found', atlasFile);
                % Needs an electrode file (even though its positions are
                % not used). 
                elecTemplateFile= fullfile(ftpath, 'template','electrode',parms.elec);
                assert(exist(elecTemplateFile,"file"), 'Electrode montage file (%s) not found', elecTemplateFile);            
            end
             
            %% Find and open the MRI from the nifti file that was created when determining
            %  the fiducials
            fprintf('***> Loading MRI data...\n');
            niftiFolder = fileparts(fids.folder);
            mriFile = fullfile(getenv('NS_ROOT'),niftiFolder,fids(1).subject + "_anat.nii");
            if exist(mriFile,"file")
                mri = ft_read_mri(char(mriFile),'dataformat','nifti');
            else
                sloc.Fiducial.dicom2nifti(fullfile(getenv('NS_ROOT'),fids.folder,fids.subject+ ".dicoms"),mriFile);
            end

            mri = ft_convert_units(mri, 'mm');

            % Align axes to CTG.(necessary for fiducial alignment below, interactive alignment
            % does this internally)
            cfg = [];
            cfg.method = 'flip';
            mri = ft_volumereslice(cfg, mri);
            cfg = [];
            cfg.method = 'fiducial';
            cfg.fiducial = fids.fiducial;
            cfg.coordsys = 'ctf'; % the desired coordinate system
            mri_realigned = ft_volumerealign(cfg, mri);
            mri_realigned.coordsys = 'ctf';

            %% Segment
            fprintf('***> Segmenting MRI data...\n');
            cfg = [];
            cfg.output    = {'gray', 'white', 'csf', 'skull', 'scalp'};
            mri_segmented  = ft_volumesegment(cfg, mri_realigned);
            mri_segmented = ft_convert_units(mri_segmented, 'mm');
            %% Head Mesh
            fprintf('***> Preparing head mesh...\n');
            cfg        = parms.mesh.cfg;
            cfg.shift  = 0;
            cfg.method = 'hexahedral';
            mesh = ft_prepare_mesh(cfg, mri_segmented);
            mesh = ft_convert_units(mesh, 'mm');

            %% Create a scalp mesh
            % Used later to project electrodes onto the scalp
            fprintf('***> Preparing scalp mesh...\n');
            cfg = parms.mesh.cfg;
            cfg.method = 'projectmesh';
            cfg.tissue = {'scalp'};
            cfg.numvertices = 10000;   % optional
            scalp = ft_prepare_mesh(cfg, mri_segmented);
            scalp = ft_convert_units(scalp, 'mm');



            %% Headmodel
            fprintf('***> Preparing head model...\n');
            cfg        = [];
            cfg.method = 'simbio';
            assert(all(string(mesh.tissuelabel)'==["csf"  "gray"  "scalp"  "skull"  "white"]),"Mesh tissue types do not match the conductivities")
            cfg.conductivity = [1.79 0.33 0.43 0.01 0.14];   % the order follows mesh.tissuelabel, which is 'csf', 'gray', 'scalp', 'skull', 'white'
            headmodel  = ft_prepare_headmodel(cfg, mesh);

            %% Sourcemodel
            fprintf('***> Preparing source model...\n');
            switch (parms.sourcemodel.cfg.method)
                case "basedonresolution"
                    % Sourcemodel based on a regular grid of dipoles
                    cfg = parms.sourcemodel.cfg;
                    
                    cfg.unit = 'mm';
                    cfg.mri = mri_segmented;
                    cfg.headmodel = headmodel;
                    cfg.headmodel.type = 'simbio';
                    sourcemodel = ft_prepare_sourcemodel(cfg);
                case "basedonmri"
                    % Sourcemodel with nodes in gray matter only.
                    cfg = parms.sourcemodel.cfg;
                    cfg.unit = 'mm';
                    cfg.mri = mri_segmented;
                    cfg.headmodel = headmodel;
                    cfg.headmodel.type = 'simbio';
                    sourcemodel = ft_prepare_sourcemodel(cfg);
                case "basedonmni"
                    % Source model based on an atlas in MNI coordinates

                    load(templateFile,'sourcemodel');
                    template = sourcemodel;clear sourcemodel;
                    template = ft_convert_units(template, 'mm');
                    thisAtlas = ft_read_atlas(char(atlasFile));
                    thisAtlas = ft_convert_units(thisAtlas, 'mm');
                    atlasVolume = thisAtlas;


                    cfg = [];
                    cfg.interpmethod = 'nearest';
                    cfg.parameter    = 'tissue';
                    thisAtlas = ft_sourceinterpolate(cfg, thisAtlas, template);

                    assert(numel(thisAtlas.tissue) == size(template.pos,1), ...
                        'Atlas volume and source-grid positions do not have matching dimensions');

                    assert(numel(template.inside) == size(template.pos,1), ...
                        'inside and source-grid positions do not have matching dimensions');

                    %% Find and snap every atlas ROI to the source grid.
                    roiIDs = find(~cellfun(@isempty, atlasVolume.tissuelabel));
                    nrROI = numel(roiIDs);
                    roiPos   = nan(nrROI,3);
                    roiLabel = cell(nrROI,1);
                    roiInside = false(nrROI,1);
                    for r = 1:nrROI
                        atlasVoxel = find(atlasVolume.tissue == roiIDs(r));
                        roiLabel{r} = atlasVolume.tissuelabel{roiIDs(r)};
                        if isempty(atlasVoxel)
                            warning('ns:Headmodel:MissingMniRoi', ...
                                'Atlas ROI %d (%s) has no voxels in the atlas volume.', ...
                                roiIDs(r), roiLabel{r});
                            continue
                        end
                        [i, j, k] = ind2sub(atlasVolume.dim, atlasVoxel);
                        atlasPos = ft_warp_apply(atlasVolume.transform, [i j k], 'homogeneous');
                        centroid = mean(atlasPos, 1);

                        % Snap the atlas centroid to the nearest actual grid point.
                        [~, nearest] = min(sum((template.pos - centroid).^2,2));
                        roiPos(r,:) = template.pos(nearest,:);
                        roiInside(r) = template.inside(nearest);
                        if ~roiInside(r)
                            warning('ns:Headmodel:MniRoiOutside', ...
                                'Snapped MNI ROI %d (%s) is outside the sourcemodel.', ...
                                roiIDs(r), roiLabel{r});
                        end
                    end
                    keep = all(isfinite(roiPos), 2);
                    roiPos = roiPos(keep,:);
                    roiLabel = roiLabel(keep);
                    roiInside = roiInside(keep);
                    nrROI = size(roiPos, 1);
                    %%
                    roiTemplate = [];
                    roiTemplate.pos      = roiPos;
                    roiTemplate.inside   = roiInside;
                    roiTemplate.outside  = ~roiInside;
                    roiTemplate.dim      = [nrROI 1 1];
                    roiTemplate.unit     = 'mm';
                    roiTemplate.coordsys  = 'mni';

                    cfg = parms.sourcemodel.cfg;
                    cfg.unit      = 'mm';                    
                    cfg.nonlinear = 'yes';
                    cfg.template  = roiTemplate;
                    cfg.mri       = mri_realigned;
                    % Elec have to be specified but their positions are not used to
                    % determine the sourcemodel. When computing the
                    % leadfield with this source model, we can specify the
                    % actual (subset of) electrodes in use in a specific
                    % experiment.
                    elec = ft_read_sens(elecTemplateFile);
                    elec = ft_convert_units(elec,'mm');
                    cfg.elec      = elec;

                    sourcemodel = ft_prepare_sourcemodel(cfg);
                    sourcemodel.tissuelabel = roiLabel;
                otherwise
                    error( 'Unknown method for preparing source model: %s.\n',parms.sourcemodel.cfg.method);
            end % switch sourcemodel.cfg.method
            fprintf('***> Uploading to the database.\n');
            tpl = key;
            tpl.mri = mri_realigned;
            tpl.headmodel = headmodel;
            tpl.scalp = scalp;
            tpl.sourcemodel = sourcemodel;
            tpl.mesh = mesh;
            tpl.segmented = mri_segmented;
            insert(tbl, tpl);
        end
    end
end
