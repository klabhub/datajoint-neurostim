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

end