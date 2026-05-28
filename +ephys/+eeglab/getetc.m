function etc = getetc(expt)
arguments
    expt 
end

if isstruct(expt)
    tbl =  ns.Experiment & expt;
elseif isa(expt,'ns.Experiment')
    tbl = expt;
end
    
cntr=0;
etc = struct;
files = (ns.File & 'extension=".mff"') &tbl;
for key  = fetch(files,'filename')'     
    mffFile = fullfile(folder(ns.Experiment &key),key.filename);
    mffFile= strrep(mffFile,'\','/'); % Avoid fprintf errors
    etcFile = strrep(mffFile,'.mff','_etc.mat');
    if exist(etcFile,'file')
        cntr= cntr+1;
        etc(cntr) = load(etcFile,'etc');
    else
        fprintf("ETC file %s does not exist.",etcFile); 
    end
end