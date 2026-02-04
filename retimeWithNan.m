function retimedT= retimeWithNan(T,newTimes,method,pv)
arguments
    T (:,:) timetable
    newTimes  %  datetimes or character vector specifying the time step
    method ="default"
    pv.TimeStep (1,1) duration = seconds(0)
    pv.SampleRate (1,1) double  = 0
    pv.EndValues (1,1) =  "extrap"
    pv.Constant (1,1)  = 123456789
    pv.IncludedEdge (1,1) string = "left"
    pv.KeepNaN (1,1) logical = true
end
keepNaN = pv.KeepNaN;
pv =rmfield(pv,'KeepNaN');
if pv.Constant == 123456789
    pv =rmfield(pv,'Constant');
end
if pv.TimeStep == seconds(0)
    pv =rmfield(pv,'TimeStep');
end
if pv.SampleRate ==0 
    pv = rmfield(pv,'SampleRate');
end
args= namedargs2cell(pv);
% First retime as instructed
retimedT= retime(T,newTimes,method,args{:});
    

if keepNaN
    % The time points that are nan in T should be nan in retimedT
    % Matlab does not have a built-in way to handle this

    % Iterate over each variable
    varNames = T.Properties.VariableNames;
    for v = 1:numel(varNames)
        col = T.(varNames{v});
        isnanCol = isnan(col);

        % Find runs of NaNs
        nanStarts = T.Time(diff([false(1,size(isnanCol,2)); isnanCol]) == 1);
        nanEnds   = T.Time(diff([isnanCol; false(1,size(isnanCol,2))]) == -1);

        %Mask bad zones in the interpolated result
        for i = 1:length(nanStarts)
            inNaNPeriod = isbetween(retimedT.Time,nanStarts(i),nanEnds(i));
            retimedT.(varNames{v})(inNaNPeriod) = NaN;
        end
    end
end

end