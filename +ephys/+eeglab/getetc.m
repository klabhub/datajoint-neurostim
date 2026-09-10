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

files = (ns.File & 'extension=".mff"') &tbl;
for key  = fetch(files,'filename')'     
    mffFile = fullfile(folder(ns.Experiment &key),key.filename);
    mffFile= strrep(mffFile,'\','/'); % Avoid fprintf errors
    etcFile = strrep(mffFile,'.mff','_etc.mat');
    if exist(etcFile,'file')
        if exist('etc','var')
            etc = catstruct(1,etc, load(etcFile,'etc').etc);
        else
            etc =   load(etcFile,'etc').etc;
        end
    else
        fprintf("ETC file %s does not exist.\n",etcFile); 
    end
end