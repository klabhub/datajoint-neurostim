function v = micronPerPixel(key)
% INPUT
%  key - A tuple containing an experiment starttime
% OUPUT
% v = micron per pixel for the x and y direction [1 2]
%

%% Get calibration info
cal = ns.getMeta(ns.Experiment &key,'xscale');
ux = str2double(unique([cal.xscale{:}]));
cal = ns.getMeta(ns.Experiment &key,'yscale');
uy = str2double(unique([cal.yscale{:}]));
assert(numel(ux)==1 && numel(uy)==1,"Magniification changed within this session??");
v = [ux uy];
end
