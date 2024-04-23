function sbx2h5(expt,pv)
% Convert sbx files to hdf5 files.
% Not yet tested (depends on sbxread from private repo).
arguments
    expt (1,1)
    pv.nrFrames (1,1) double  = Inf
    pv.blockSize (1,1) double =200
end

if isa(expt,'ns.Experiment')
    expt = fetch(expt);
end

info  =sbx.readInfoFile(expt);
sz = info.sz([2 1]); % Swap x/y dims to mathc permute below

sbxFile = ((ns.File & expt) & 'extension=''.sbx''');
fname = fetch1(sbxFile & expt,'filename');
fldr = folder(ns.Experiment & expt);
ff =fldr + fname(1:end-3) + 'mat';
if ~exist(ff,'file')
    error('Sbx File %s not found',strrep(ff,'\','/'))
end

h5File = fldr + fname(1:end-3)  + ".h5";
z = sbxread(sbxFile,1,1);


pv.nrFrames = min(pv.nrFrames,info.max_idx);

framesReadSoFar=0;
framesToRead = min(pv.blockSize,pv.nrFrames-framesReadSoFar);
while(framesToRead>0)
    try
        q = sbxread(sbxFile,framesReadSoFar,framesToRead);
        q = squeeze(q(1,:,:,:)); % extract green channel only
        q = permute(q,[2 1 3]);
        if(framesReadSoFar==0)
            % First block; create the file
            h5create(h5File,'/data',[sz Inf],'DataType','uint16','ChunkSize',[sz framesToRead]);
        end
        % Add to the h5 file
        h5write(h5File,'/data',q,[1 1 framesReadSoFar+1],[sz framesToRead]);
    catch
        framesToRead=0;
    end
    framesReadSoFar = framesReadSoFar+framesToRead;
    framesToRead = min(pv.blockSize,pv.nrFrames-framesReadSoFar);
end

