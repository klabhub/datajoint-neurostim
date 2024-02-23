function key = stripToPrimary(tbl,key)
% Based on a key struct that may have additional fields, strip all fields that are
% not part of this tables primary key. 
% INPUT
% tbl  - Relvar table
% key - Struct with (a subset of) primary keys as values (plus any  additional)
% fields).
% OUTPUT
% pk - Struct with only keys that are primary key values.
%
pk  = tbl.primaryKey;
notPK  =setdiff(fieldnames(key),pk);
if ~isempty(notPK)
    key = rmfield(key,notPK);
end
end
            