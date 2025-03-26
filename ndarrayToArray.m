function m = ndarrayToArray(x, pv)
% Convert a numpy array to a matlab array, with quick conversions for
% matrices, and slower conversions for ragged arrays/
arguments
    x (1,1) py.object
    pv.single (1,1) logical = false
end

dtype = string(x.dtype.name);
numericTypes = ["float64", "float32", "int64", "int32", "int16", "uint8", "bool"];
isNumericArray = ismember(dtype, numericTypes) && ~strcmp(dtype, "object");

if isNumericArray
    % Get shape
    shape = double(x.shape);

    % Flatten the numpy array and convert to Python array.array
    if pv.single
        pyarr = py.array.array('f', x.astype('float32').flatten().tolist());
        m = single(pyarr);
    else
        pyarr = py.array.array('d', x.astype('float64').flatten().tolist());
        m = double(pyarr);
    end

    % Reshape back to original shape
    m = reshape(m, fliplr(shape))';
    return
end

% Fallback for object/ragged arrays
asCell = cell(x.tolist());
op = ternary(pv.single, @single, @double);
rowLengths = cellfun(@numel, asCell);
isUniform = all(rowLengths == rowLengths(1));

if isUniform
    m = cellfun(op, asCell, 'UniformOutput', false);
    m = cell2mat(m');
else
    warning("Ragged structure detected; returning cell array.");
    m = cellfun(op, asCell, 'UniformOutput', false);
end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
