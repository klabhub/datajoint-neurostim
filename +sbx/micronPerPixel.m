function v = micronPerPixel(key)
% INPUT
%  key - A tuple containing an experiment starttime
% OUPUT
% v = micron per pixel for the x and y direction [1 2]
%

%% Get calibration info
cal = fetch(ns.ExperimentMeta &key & struct('meta_name','xscale'),'meta_value');
ux = unique([cal.meta_value]);
cal = fetch(ns.ExperimentMeta &key & struct('meta_name','yscale'),'meta_value');
uy = unique([cal.meta_value]);
assert(numel(ux)==1 && numel(uy)==1,"Magniification changed within this session??");
v = [ux uy];
end
