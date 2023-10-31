function out = getMeta(tbl,meta,pv)
% Convenience function to retrieve meta data for Subject, Session, or
% Experiment tables.
% 
% INPUT
% tbl  -  A subject,session, or experiment table
% meta - name(s) of the meta parameter to retrieve
% OUTPUT
% T  -  A Matlab Table with columns for the primary keys of the table, plus columns
%       for each of the requested meta parameters. Note that meta data are
%       strings.
arguments
    tbl (1,1) {mustBeA(tbl,{'ns.Experiment','ns.Subject','ns.Session'})}
    meta {mustBeText} = ""
    pv.type (1,1) string = "table"
    pv.number = inf
end
% Get the meta table
metaTable = feval([class(tbl) 'Meta']);

if strlength(meta)==0
    % get all meta fields for this tbl
    meta = unique({fetch(metaTable & tbl,'meta_name').meta_name});
elseif ischar(meta)
    % Make cellstring for loop
    meta = {meta};
end
% Loop over meta fields to add
T = fetchtable(tbl);
for i=1:numel(meta)
    thisMetaT= fetchtable(metaTable &tbl & struct('meta_name',meta{i}),'meta_value');
    thisMetaT = addvars(thisMetaT,thisMetaT.meta_value,'NewVariableNames',meta{i});  
    if ~isempty(thisMetaT)
        T = outerjoin(T,thisMetaT,'MergeKeys',true,'RightVariables',setdiff(thisMetaT.Properties.VariableNames,{'meta_name','meta_value'}));    
    end
end


switch (pv.type)
    case 'table'
        % nothing to do
        out = T;
    case 'double'
        % Meta is a number stored as string.
        out = cellfun(@isempty,T{:,meta});
        if any(~out)
            out =cellfun(@str2num,T{~out,meta});
        else
            error('No %s found in this tbl',meta{:})
        end
    otherwise 

end




