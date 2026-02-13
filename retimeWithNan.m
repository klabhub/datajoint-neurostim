function retimedT= retimeWithNan(T,newTimes,method,pv)
% This function wraps around the built-in retime, with an additional option
% (AnyNanIsNanSum) that specifies how to handle missing data with the "sum"
% retime method
% If there are no missing data, this function just calls retime
% If there are missing data:
%   AnyNanIsNanSum=true, a bin with one ore more NaN and method="sum" will  result in NaN.  
%   AnyNanIsNanSum=false, a bin with all NaN result in NaN, but all other bins will result in the sum of the non-NaN elements
%
% The default is  AnyNanIsNanSum =true.
arguments
    T (:,:) timetable
    newTimes  %  datetimes or character vector specifying the time step
    method ="linear"
    pv.TimeStep (1,1) duration = seconds(0)
    pv.SampleRate (1,1) double  = 0
    pv.EndValues (1,1) =  "extrap"
    pv.Constant (1,1)  = 123456789
    pv.IncludedEdge (1,1) string = "left"
    pv.AnyNanIsNanSum (1,1) logical = true
end

% Massage the pv so that we can pass this to the builtin retime
anyNanIsNaNSum = pv.AnyNanIsNanSum;
pv =rmfield(pv,'AnyNanIsNanSum');
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

if method=="sum" 
    % How to handle NaN?    
    % Matlab's default for sum is to simply ignore NaN, which results in 0 for a
    % time bin with all NaN.       
    if anyNanIsNaNSum    
        % Even just one nan in the bin results in a NaN
        method = @(x) sum(x, 'includemissing'); 
    else
        % Only all-Nan results in NaN
        % For a timetable without NaN, this is the same as the built-in
        % retime
        method = @(x) sum(x, 'omitmissing') + 0./~all(isnan(x)); 
    end    
end
% Call the built-in retime
retimedT= retime(T,newTimes,method,args{:});
  
end