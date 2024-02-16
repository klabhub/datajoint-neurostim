function  [out]  = compare(tbl, r1,r2,pv)
arguments
    tbl (1,1)   % A dj.RelVar (table)
    r1  struct % a DJ restriction 
    r2  struct % another DJ restriction    
    pv.variables = tbl.nonKeyFields;    
    pv.keys = tbl.primaryKey;
end
    if ischar(pv.variables)
        pv.variables = {pv.variables};
    end
t1 =fetchtable(tbl&r1,pv.variables{:});
t2 =fetchtable(tbl&r2,pv.variables{:});

r1Names = fieldnames(r1);
r2Names = fieldnames(r2);

T= innerjoin(t1,t2,"Keys",setdiff(pv.keys,union(r1Names,r2Names)));


if nargout==0
    tiledlayout("flow")
    for i=1:numel(pv.variables)
        nexttile
        x = T.(pv.variables(i) + "_t1");
        y =T.(pv.variables(i) + "_t2");
        plot(x,y,'.')
        xlabel (pv.variables(i) + " " + strjoin(namedargs2cell(r1),"/"))
        ylabel (pv.variables(i) + " " + strjoin(namedargs2cell(r2),"/"))
        lims = [min([x;y]) max([x;y])];
    xlim(lims)
    ylim(lims)
    hold on
    plot(xlim,ylim,'k')
    end
else
    out = T;
end