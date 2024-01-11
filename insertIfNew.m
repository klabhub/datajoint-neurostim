function insertIfNew(tbl,tpl)
% Check whether a tpl already exists in the table and insert it if not. 

arguments
    tbl (1,1) 
    tpl (:,:) struct
end

for i=1:numel(tpl)
    key = ns.stripToPrimary(tbl,tpl(i));
    if count(tbl&key)==0
        insert(tbl,tpl(i))
    end
end
