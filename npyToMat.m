function npyToMat(npyFile,matFile)
% Function that saves an .npy file as a .mat file, using scipy.io.savemat.
% The actual code runs in Python, using pyrun.
% INPUT
% npyFile   - the source .npy file
% matFile  - the target .mat file (defaults to the same file name as the
%               .npy file but with a .mat extension).
% BK - Mar 2023.
arguments
    npyFile {mustBeText}
    matFile {mustBeText} = ""
end
if matFile==""
    [fldr,file,~] = fileparts(npyFile);
    matFile = fullfile(fldr,string(file) + ".mat");
end
pyCode = 'from scipy.io import savemat; import numpy as np;stat = np.load(''%s'',allow_pickle=True);savemat(''%s'',mdict={''stat'': stat})';
pyCode = sprintf(pyCode,strrep(npyFile,'\','/'), strrep(matFile,'\','/'));
pyrun(pyCode);
end
