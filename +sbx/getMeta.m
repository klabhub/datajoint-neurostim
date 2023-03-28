function T = getMeta(tbl,meta)
% Convenience function to retrieve meta data for Subject, Session, or
% Experiment tables.
% 
% INPUT
% tbl  -  A subject,session, or experiment table
% meta - name(s) of the meta parameter to retrieve
% OUTPUT
% T  -  A Matlab Table with columns for the primary keys of the table, plus columns
%       for each of the requested meta parameters. 
arguments
    tbl (1,1) {mustBeA(tbl,{'ns.Experiment','ns.Subejct','ns.Session'})}
    meta {mustBeText}
end
if ischar(meta)
    meta = {meta};
end
% Consider only files that don't have the info meta data 
metaTable = feval([class(tbl) 'Meta']);
T = fetchtable(tbl & metaTable & struct('meta_name',meta));

for i=1:numel(meta)
    thisMetaT= fetchtable(metaTable &tbl & struct('meta_name',meta{i}),'meta_value');
    thisMetaT = addvars(thisMetaT,thisMetaT.meta_value,'NewVariableNames',meta{i});
    thisMetaT =removevars(thisMetaT,["meta_name","meta_value"]);
    T = innerjoin(T,thisMetaT);
end

