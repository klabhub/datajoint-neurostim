function chunkedInsert(tbl,tpl)
% Insert the struct array of tpls into datajoint table tbl.
% The number of tpls is chosen to send approximately NS_BYTESPERINSERT
% at a time. This value is read from the environment
% (getenv("NS_BYTESPERINSERT")) or defaults to 5e6 bytes.
bytesPerInsert = getenv("NS_BYTESPERINSERT");
if isempty(bytesPerInsert)
    bytesPerInsert = 50e6; % 50MB default chunk - ok for DJ server with large innodb_buffer and ram.
end

get_mem_size = @(x) whos('x').bytes;

% Determine maximum chunk size based on individual element size
totalMem = get_mem_size(tpl);
elementSize = totalMem / numel(tpl);
maxElementsPerChunk = floor(bytesPerInsert / elementSize);

% Handle cases where individual element is too large.
if maxElementsPerChunk < 1
    maxElementsPerChunk = 1;  % Ensure at least one element is sent.
    fprintf('Individual struct element exceeds bytesPerInsert. Sending one element at a time.\n');
end


tic;
nrTpls = numel(tpl);


fprintf('Uploading to server (%d tuples with %d MB per chunk):',maxElementsPerChunk,round(maxElementsPerChunk*elementSize/1e6))

i = 1;
progress_old = 0;

while i <= nrTpls

    % Find optimal chunk size for this iteration without exceeding limit
    currentChunkSize = min(maxElementsPerChunk, nrTpls - i + 1);    
    thisChunk = i:min(nrTpls,(i+currentChunkSize-1));
    insert(tbl,tpl(thisChunk));
    i = i + currentChunkSize;  % Move to the next chunk
    progress_new = round(100*(i-1)/nrTpls);
    if (progress_new - progress_old) >= 10

        fprintf("\t %d%%", progress_new);
        progress_old = progress_new;
    end
end

fprintf('\n\tDone in %d seconds. Uploaded %d MB \n.',round(toc),round(totalMem/1e6));
end