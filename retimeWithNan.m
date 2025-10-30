function retimedT= retimeWithNan(T,newTimes,pv)
arguments
    T (:,:) timetable
    newTimes (:,1) duration
    pv.interpolation (1,1) string = "linear"
    pv.endValues (1,1) =  NaN
    pv.keepNaN (1,1) logical = true
end

% First retime as instructed
retimedT= retime(T,newTimes,pv.interpolation,'EndValues',pv.endValues);

if pv.keepNaN
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