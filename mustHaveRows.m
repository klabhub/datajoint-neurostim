function mustHaveRows(tbl,nrRequired)
%MUSTHAVEROWS Validate that a DataJoint relation contains required rows.
%   With no required count, the relation must be nonempty. With a required
%   count, the relation must contain exactly that many rows. Counting is
%   bounded at the requested count, so validation does not scan a large
%   relation just to reject it.

if nargin < 2
    if ~hasAtLeast(tbl,1)
        throwAsCaller(MException('mustHaveRows:zeroRows', ...
            'This DataJoint table has no rows'))
    end
    return
end

[nrRows,isCapped] = countUpTo(tbl,nrRequired);
if isCapped || nrRows ~= nrRequired
    if isCapped
        actualRows = sprintf('more than %d',nrRequired);
    else
        actualRows = sprintf('%d',nrRows);
    end
    throwAsCaller(MException('mustHaveRows:mismatchedRows', ...
        'This DataJoint table has %s rows (%d required)', ...
        actualRows,nrRequired))
end
end
