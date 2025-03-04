function m = ndarrayToArray(x,pv)
arguments
    x (1,1) py.numpy.ndarray 
    pv.single (1,1) logical  = false
end

m = cell(x.tolist());
if pv.single
op = @single;
else
    op = @double;
end
m = cellfun(op,m,UniformOutput=false);
m = cell2mat(m');
                 