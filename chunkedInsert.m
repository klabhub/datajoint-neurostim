function chunkedInsert(tbl,tpl)
% Insert the struct array of tpls into datajoint table tbl.
% The number of tpls is chosen to send approximately NS_BYTESPERINSERT
% at a time. This value is read from the environment
% (getenv("NS_BYTESPERINSERT")) or defaults to 5e6 bytes.
bytesPerInsert = getenv("NS_BYTESPERINSERT");
if isempty(bytesPerInsert)
    bytesPerInsert = 5e6; % 5MB default chunk
end

get_mem_size = @(x) whos('x').bytes;

% Determine maximum chunk size based on individual element size
totalMem = get_mem_size(tpl);
elementSize = totalMem / numel(tpl);
maxElementsPerChunk = floor(bytesPerInsert / elementSize);

% Handle cases where individual element is too large.
if maxElementsPerChunk < 1
    maxElementsPerChunk = 1;  % Ensure at least one element is sent.
    warning('Individual struct element exceeds bytesPerInsert. Sending one element at a time.');
end


tic;
nrTpls = numel(tpl);
fprintf('Uploading to server: ')

i = 1;
progress_old = 0;

while i <= nrTpls

    % Find optimal chunk size for this iteration without exceeding limit
    currentChunkSize = min(maxElementsPerChunk, nrTpls - i + 1);

    % Double check if the submission exceeds the limit or not
    thisChunk = i:(i+currentChunkSize-1);
    while (get_mem_size(tpl(thisChunk)) > bytesPerInsert) && currentChunkSize >1

        currentChunkSize = currentChunkSize - 1;
        thisChunk = i:(i+currentChunkSize-1);

    end

    insert(tbl,tpl(thisChunk));
    
    i = i + currentChunkSize;  % Move to the next chunk
    progress_new = round(100*(i-1)/nrTpls);
    if (progress_new - progress_old) >= 10

        fprintf("\t %d%%", progress_new);
        progress_old = progress_new;
    end
end

fprintf('\n\tDone in %d seconds.\n.',round(toc))
end