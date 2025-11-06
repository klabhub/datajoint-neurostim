function tpl =makeMymSafe(tpl)
% Take an array of tuples and recursively convert strings to chars and
% function_handles to function names. This is necessary to store the tuple
% in a datajoint database using MyM.
arguments
    tpl (:,1) struct
end

for t = 1:numel(tpl)
    fn =fieldnames(tpl(t));
    nrFields= numel(fn);
    for f=1:nrFields
        thisField= tpl(t).(fn{f});
        if isobject(thisField) && ~isstring(thisField)
            warning('off','MATLAB:structOnObject')
            thisField = struct(thisField);
            warning('on','MATLAB:structOnObject')
            tpl(t).(fn{f}) = thisField;
        end
        switch class(thisField)
            case 'struct'
                    tpl(t).(fn{f}) = makeMymSafe(thisField); % Recurse
            case 'cell'
                    isStr  = cellfun(@isstring,thisField);
                    tpl(t).(fn{f})(isStr) =cellfun(@char,thisField(isStr),'UniformOutput',false);
                    isStruct =cellfun(@isstruct,thisField);
                    tpl(t).(fn{f})(isStruct) =cellfun(@makeMymSafe,thisField(isStruct),'UniformOutput',false);
            case 'string'
                     tpl(t).(fn{f})  = char( thisField);
            case 'function_handle'
                    tpl(t).(fn{f}) = func2str(thisField);
            otherwise
                %Nothing to do
        end
    end  
end
end