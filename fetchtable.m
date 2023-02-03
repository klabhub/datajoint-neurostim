function T = fetchtable(varargin)
% Wrapper around the fetch function in DataJoint to return a Matlab table
% instead of a struct array.

tuples= fetch(varargin{:});
try
T = struct2table(tuples);
catch
    T=tuples;
end