function T = fetchtable(varargin)
% Wrapper around the fetch function in DataJoint to return a Matlab table
% instead of a struct array.

tuples= fetch(varargin{:});
T = struct2table(tuples,'AsArray',true);
end