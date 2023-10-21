function v = micronPerPixel(key)
% INPUT
%  key - A tuple containing an experiment starttime
% OUPUT
% v = micron per pixel for the x and y direction [1 2]
%

%% Get calibration info
ux  = unique(ns.getMeta(ns.Experiment &key,'xscale',type="double"));
uy = unique(ns.getMeta(ns.Experiment &key,'yscale',type="double"));
assert(numel(ux)==1 && numel(uy)==1,"Magniification changed within this session??");
v = [ux uy];
end
