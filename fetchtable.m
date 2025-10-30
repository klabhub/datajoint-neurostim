function T = fetchtable(tbl,args)
% Wrapper around the fetch function in DataJoint to return a Matlab table
% instead of a struct array.
arguments
    tbl 
end
arguments (Repeating)
    args (1,:) {mustBeText} % Arguments to pass to fetch
end
% Convert args to a cell array of char vectors
if nargin == 2 && (isstring(args{1}) && numel(args{1}) > 1)
    % Single string vector argument
    charArgs = cellstr(args{1});
else
    % Multiple arguments or single string/char
    charArgs = cellfun(@char, args, 'UniformOutput', false);
end

% Call fetch with table and converted arguments
tuples= fetch(tbl,charArgs{:});
T = struct2table(tuples,'AsArray',true);
T= convertvars(T,@iscellstr,"string");
end