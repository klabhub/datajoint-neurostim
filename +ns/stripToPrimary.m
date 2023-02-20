function stripped = stripToPrimary(tbl,key)
% Based on a key struct that may have additional fields, strip all fields that are
% not part of this tables primary key. 
% INPUT
% tbl  - Relvar table
% key - Struct with (a subset of) primary keys as values (plus any  additional)
% fields).
% OUTPUT
% stripped - Struct with only keys that are primary key values.
%
pkNames  = tbl.primaryKey;
out = setdiff(fieldnames(key),pkNames);
stripped = rmfield(key,out);
            