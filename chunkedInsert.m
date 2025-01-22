function chunkedInsert(tbl,tpl)
% Insert the struct array of tpls into datajoint table tbl.
% The number of tpls is chosen to send approximately NS_BYTESPERINSERT
% at a time. This value is read from the environment
% (getenv("NS_BYTESPERINSERT")) or defaults to 5e6 B.


bytesPerInsert = getenv("NS_BYTESPERINSERT");
if isempty(bytesPerInsert)
    bytesPerInsert = 5e6; % 5MB default chunk
end
totalMem =whos("tpl").bytes;
chunkSize = floor(totalMem/bytesPerInsert); 
tic;
nrTpls = numel(tpl);
fprintf('Uploading to server ')
for i=1:chunkSize:nrTpls
    fprintf('.')
    if mod(i,80)==0;fprintf('\n');end
    thisChunk = i:min(nrTpls,i+chunkSize-1);
    insert(tbl,tpl(thisChunk));
end
fprintf('Done in %d seconds.\n.',round(toc))
end