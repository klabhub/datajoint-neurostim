function qry = in(column, list)
% Convenience function to add an IN where SQL query
% For instance (in("type" ,["A","B","C"]), will generate the query string
% "type IN('A','B','C')"
% which can be added as a string to a datajoint query using &.
arguments
    column (1,1) string
    list (1,:) string
end
qry  = column + " IN('" + strjoin(list,''',''') + "')";
end