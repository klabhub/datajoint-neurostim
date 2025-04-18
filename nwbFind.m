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
% to the ns package. Force loading all requested packages (TODO : figure out how to
% determine all packages in the database schema and load all of those by default)
ret = dj.ERD();
if packages~=""
    for package = packages
        obj = dj.ERD(feval(package + ".getSchema"));
        ret = ret + obj;
    end
end

className = string([]);
try
cls = feval(top);
m = metaclass(cls);
if ismember('nwb',{m.MethodList.Name})
    % NWB export defined
    className = m.Name;
end

children=cls.children;
for ch=1:numel(children)
    childName =  string(dj.conn().tableToClass(children{ch}));
    className = [className nwbFind(childName)]; %#ok<AGROW>
end
catch me
    me.message

end