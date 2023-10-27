function insertIfNew(tbl,tpl)
% Check whether a tpl already exists in the table and insert it if not. 

arguments
    tbl (1,1) 
    tpl (1,1) struct
end

key = ns.stripToPrimary(tbl,tpl);
if count(tbl&key)==0
    insert(tbl,tpl)
end
