function out = readInfoFile(expt)
% Given an expt tuple, read and return the info struct from the sbx .mat output file.
if isa(expt,'ns.Experiment')
    expt = fetch(expt);
end

% Loop over experiments
for e=1:numel(expt)
    sbxFile = ((ns.File & expt(e)) & 'extension=''.sbx''');
    % There should be exactly 1 sbx file per expt. We had 3 in at least one instance
    nrSbxFiles = count(sbxFile);
    if nrSbxFiles ~=1
        error("%d sbx files in experiment (%s on %s @ %s)",count(sbxFile),expt(e).subject,expt(e).session_date,expt(e).starttime);
    end

    fname = fetch1(sbxFile & expt(e),'filename');
    fldr = folder(ns.Experiment & expt(e));
    ff =fldr + fname(1:end-3) + 'mat';
    if exist(ff,'file')
        load(ff ,'info'); % Get the info struct
    else
        error('Could not load sbx info struct. File %s not found',strrep(ff,'\','/'))
    end


    switch info.scanbox_version
        case 3
            info.xscale = info.dxcal;
            info.yscale = info.dycal;

            % % Calibration contains curvefit objects, which mym does not like.
            % Remove/replace
            if exist('fittype',"class")
                % If @fittype exists, we can use curve fit toolbox to
                % extract the functions as strings and parameters as values
                info.calibration.formulaX = formula(info.calibration.fx);
                info.calibration.formulaY = formula(info.calibration.fy);
                info.calibration.parmsX = [info.calibration.fx.a info.calibration.fx.b info.calibration.fx.c info.calibration.fx.d];
                info.calibration.parmsY = [info.calibration.fy.a info.calibration.fy.b info.calibration.fy.c info.calibration.fy.d];
            else
                fprintf(2,"The curve fitting toolbox is needed for calibration. Trying to extract from struct");
                warning( 'off','MATLAB:structOnObject')
                tmp = struct(info.calibration.fx);
                info.calibration.formulaX = tmp.defn;
                info.calibration.parmsX = [tmp.coeffValues{:}];
                tmp = struct(info.calibration.fy);
                info.calibration.formulaY = tmp.defn;
                info.calibration.parmsY = [tmp.coeffValues{:}];     
                warning( 'on','MATLAB:structOnObject')
            end
            info.calibration = rmfield(info.calibration,'fx');
            info.calibration = rmfield(info.calibration,'fy');        %            
        otherwise
            cal = info.calibration(info.config.magnification);
            info.xscale = 1./cal.x;
            info.yscale = 1./cal.y;
    end

    % Determine the number of planes  from the optotune settings.
    if isfield(info,'otwave')
        nrPlanes = max(1,numel(info.otwave));
    elseif isnan(info.otparam(3))
        nrPlanes = 1 ;
    else
        nrPlanes = info.otparam(3);
    end
    nrChannels = (info.channels==0)+1;  % 0 means red+green, all others mean single channel
    % 2 bytes per pixel.
    dirInfo = dir(fullfile(fldr,fname));
    info.nrFrames = dirInfo.bytes./prod(info.sz)/nrChannels/nrPlanes/2;
    info.nrPlanes  = nrPlanes;
    if floor(info.nrFrames) ~=info.nrFrames
        % This seems to happen in 2-plane recordings;  an extra byte in the
        % file?
        fprintf(2,"Half frames (nrPlanes = %d) in %s?? Flooring.\n",nrPlanes,ff)
        info.nrFrames = floor(info.nrFrames);
    end  
    out(e) = info; %#ok<AGROW>
end
end