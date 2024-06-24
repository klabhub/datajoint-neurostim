function className = nwbFind(top,packages)
% Starting from the top of a schema, walk down the tree
% to identify all tables that have a method called  nwb.
% This method is reserved for Neurodata Without Borders export.
% See ns.Experiment/nwbExport
arguments
    top (1,1) string = "ns.Subject"
    packages (1,:) string = ""
end

% Unless a package has already been used in a session, children are limited
% to the ns package. Force loading all packages (TODO : figure out how to
% determine all packages in the database schema)
if packages ~=""
    currentPackages=dj.conn().packages.keys;
    for p=packages
        if ismember(p,currentPackages);continu;end
        dbase = dj.conn().schemas.keys;
        dj.Schema(dj.conn,p,dbase{1});
    end
end
className = string([]);
cls = feval(top);
m = metaclass(cls);
if ismember('nwb',{m.MethodList.Name})
    % NWB export defined
    className = m.Name;
end

children=cls.children;
for ch=1:numel(children)
    childName =  string(dj.conn().tableToClass(children{ch}));
    className = [className findNwb(childName)]; %#ok<AGROW>
end