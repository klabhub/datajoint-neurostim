function classNames = nwbFind(top,packages)
% Starting from the top of a schema, walk down the tree
% to identify all tables that have a method called  nwb.
% This method is reserved for Neurodata Without Borders export.
% See ns.Experiment/nwbExport
% 
% Class names will be sorted by the order of the packages and then
% alphabetically within a package.
% In other words, if packages = ["ns" "sbx"] , all ns. Classes will be
% listed first.  Ordering is required because some elements of the NWB
% hierarchy are setup in one place, then used in another.
arguments
    top (1,1) string = "ns.Subject"
    packages (1,:) string = ""
end

classNames = string([]);

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

try
    cls = feval(top);
    m = metaclass(cls);
    if ismember('nwb',{m.MethodList.Name})
        % NWB export defined
        classNames = m.Name;
    end

    children=cls.children;
    for ch=1:numel(children)
        childName =  string(dj.conn().tableToClass(children{ch}));
        classNames = [classNames nwbFind(childName)]; %#ok<AGROW>
    end
catch me
    me.message
end

if packages ~=""
    % Sort such that the members of the first package occur first.
    sortedClassNames = [];
    for i=1:numel(packages)
        sortedClassNames = [sortedClassNames, sort(classNames(startsWith(classNames,packages(i) + ".")))]; %#ok<AGROW>
    end
    classNames = sortedClassNames;
end
end