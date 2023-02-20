function v = datetime(tpl)
% Extract a datetime from a tpl that contains session_date and
% starttime
v = unique(datetime(string({tpl.session_date}) + " " + string({tpl.starttime})));
end