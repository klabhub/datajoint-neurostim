function info = readInfoFile(expt)
% Given an expt tuple, read and return the info struct from the sbx .mat output file.

sbxFiles = ((ns.File & expt) & 'extension=''.sbx''');
fname = fetch1(sbxFiles & expt,'filename');
fldr = folder(ns.Experiment & expt);
ff =fldr + fname(1:end-3) + 'mat';
if exist(ff,'file')
    load(ff ,'info'); % Get the info struct
else
    error('Could not load sbx info struct. File %s not found',strrep(ff,'\','/'))
end


if info.scanbox_version==3
    % % Calibration contains curvefit objects, which mym does not like.
    % Remove/replace
    info.calibration.formulaX = formula(info.calibration.fx);
    info.calibration.formulaY = formula(info.calibration.fy);
    info.calibration.parmsX = [info.calibration.fx.a info.calibration.fx.b info.calibration.fx.c info.calibration.fx.d];
    info.calibration.parmsY = [info.calibration.fy.a info.calibration.fy.b info.calibration.fy.c info.calibration.fy.d];
    info.calibration = rmfield(info.calibration,'fx');
    info.calibration = rmfield(info.calibration,'fy');
end

end