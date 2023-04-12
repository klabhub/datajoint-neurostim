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


% Create a single cell with parameter value pairs.


pvPairs = cell(1,numel(fn)*2);
[pvPairs{1:2:end}] =deal(fn{:});
[pvPairs{2:2:end}] =deal(vals{:});
% Create the merged struct
v =struct(pvPairs{:});
end