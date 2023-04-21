function x1 = catpad(x1,x2)
% Add a vector (x2) as a column in a matrix (x1). x1 and x2 are allowed to have
% different number of rows, missing rows will be padded with nans at the end.
% Hence the assumption is that the start point (row 1) for x1 and x2
% corresponds to the same sample, but x2 can have more or fewer samples at
% the end.
%
% INPUT
% x1 - Matrix/Vector 
% x2 - Matrix/Vector

if isempty(x1)
    x1 =x2; 
else
    n1 = size(x1,1);
    n2 = size(x2,1);
    pad1 = max(n2-n1,0);
    pad2 = max(n1-n2,0);
    sz1 = size(x1);
    sz1(1) = pad1;
    sz2 = size(x2);
    sz2(1) = pad2;
    x1 = [[x1;nan(sz1)] [x2;nan(sz2)]];
end
end