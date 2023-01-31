function thisJson= readJson(file)
% Utility to read json files and change all numeric values to
% strings
fid =fopen(file,'r');
txt = fread(fid,inf,'*char');
fclose(fid);
thisJson = jsondecode(txt');
fn = fieldnames(thisJson);
for j=1:numel(thisJson)
    for i=1:numel(fn)
        if isnumeric(thisJson(j).(fn{i})) || ischar(thisJson(j).(fn{i}))
            thisJson(j).(fn{i}) = string(num2str(thisJson(j).(fn{i})));
        end
    end
end
end
