function qry = in(column, list)
arguments
    column (1,1) string
    list (1,:) string
end

qry  = column + " IN('" + strjoin(list,''',''') + "')";
end