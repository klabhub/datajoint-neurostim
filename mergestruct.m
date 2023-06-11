function v= mergestruct(varargin)
% Create a struct by combining the fields in all input structs.
% If one of the input structs is a singleton and the other(s) an array, the
% scalar struct will be repmatted (i.e. used in the entire array).
%
% INPUT
% varargin = Comma separated list of structs.
% OUTPUT
% v = The merged struct
%
% BK - Dec 2022

nrStructs = numel(varargin);
arraySize= cellfun(@numel,varargin);
singleton = arraySize==1;
expandTo  = unique(arraySize(~singleton));
assert(numel(expandTo)==1,'mergestruct can only merge singletons with one struct array size');
[varargin(singleton)] =cellfun(@(x) repmat(x,[expandTo 1]),varargin(singleton),'uni',false);

% Extract fieldnames and values from each input struct
fn = {};
vals ={};
for i=1:nrStructs
    newFields = fieldnames(varargin{i});
    newVals = struct2cell(varargin{i});
    % Check for duplicate
    [duplicate,ix] = ismember(newFields, fn);
    if any(duplicate)
        same = cellfun(@isequaln,newVals(duplicate),vals(ix(duplicate)));
        if  all(same)
            % Duplicates with identical values- ignore
            newFields = newFields(~duplicate);
            newVals = newVals(~duplicate);
        else
            nf = newFields(duplicate);
            nv = newVals(duplicate);
            ov = vals(ix(duplicate));
            field = nf(~same)
            newValue  = nv(~same) 
            oldValue = ov(~same)
            error('Duplicate fields with mismatched values')
        end
    end

    fn =cat(1,fn,newFields);
    vals = cat(1,vals,newVals);
end


% Create the merged struct.
v= cell2struct(vals,fn);
end