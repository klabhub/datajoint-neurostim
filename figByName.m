function out = figByName(name,props)
% Find a a figure by name, bring it to the front , or create it if it does
% not yet exist. An optional struct with figure properties can be passed. 
arguments
    name (1,1) string 
    props (1,1) struct = struct;
end
fig = findobj(get(0,'children'),'flat','name',name);
if isempty(fig)
    fig = figure;
    set(fig,'Name',name);
else
    figure(fig);
end

if ~isempty(props)
    set(fig,props)
end

if nargout>0
    out =fig;
end