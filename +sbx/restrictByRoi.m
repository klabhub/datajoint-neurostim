function r = restrictByRoi(tbl,varargin)

sessions  = ns.Session & tbl;
allTpls=  [];
for sess = fetch(sessions)'
    thisTbl = sbx.PreprocessedRoi &sess;
    thisTbl.restrict(varargin{:});
    keep= fetch(thisTbl,'roi');
    thisRois = num2cell([keep.roi]);
    thisTpls = fetch((tbl & sess) & struct('channel',thisRois'));
    if isempty(allTpls)
        allTpls =thisTpls;
    else
        allTpls = cat(1,allTpls,thisTpls);
    end
end

r =tbl & allTpls;