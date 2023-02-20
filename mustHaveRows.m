function mustHaveRows(tbl,nrRequired)
% Validation function to check the number of required rows in a DataJoint
% table.  
% INPUT
% tbl  - The DJ table
% nrRequired - The number of rows that is required . Defaults to >0
% OUTPUT
% void - To be used in a function arguments validation block:
% arguments
% tbl {mustHaveRows(tbl,1)}  - To allow only a single row table.

nrRows = count(tbl);
if nargin <2 
    if nrRows==0    
        throwAsCaller(MException('mustHaveRows:zeroRows','This DataJoint table has no rows'))
    end
else
    if nrRows~=nrRequired
        throwAsCaller(MException('mustHaveRows:mismatchedRows',sprintf('This DataJoint table has %d rows (%d required)',nrRows,nrRequired)))
    end
end
end