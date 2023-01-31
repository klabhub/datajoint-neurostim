function v= mergestruct(varargin)
% Create a scalar struct by combining the fields in all input structs.
%
% INPUT
% varargin = Comma separated list of structs.
% OUTPUT
% v = A struct
%
% BK - Dec 2022

nrStructs = numel(varargin);

% Extract fieldnames and values from each input struct
fn = {};
vals ={};
for i=1:nrStructs
    fn =cat(1,fn,fieldnames(varargin{i}));
    vals = cat(1,vals,struct2cell(varargin{i}));
end
% Create a single cell with parameter value pairs.

pvPairs = cell(1,numel(fn)*2);
[pvPairs{1:2:end}] =deal(fn{:});
[pvPairs{2:2:end}] =deal(vals{:});
% Create the merged struct
v =struct(pvPairs{:});
end