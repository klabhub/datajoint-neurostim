function [v,t] = timetableToDouble(T)
% Converts a timetable T to a 3D array of signal values and the
% corresponding time points.
arguments
    T timetable
end

[nrTimePoints, nrTrials] = size(T);
if nrTimePoints>0
    nrChannels = numel(T{1,1});
    v = permute(reshape(T.Variables,[nrTimePoints nrChannels nrTrials]),[1 3 2]);
    t = T.Time;
else
    v=[];
    t=[];
end