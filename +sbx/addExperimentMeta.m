function addExperimentMeta(key,pv)
% After running nsScan to populate the ns.Experiments and ns.File table,
% this funciton will add the metadata from the scanbox .mat file (i.e. the 
% the number of frames (nrframes)  and the scaling of x and y pixels ) to the
% ExperimentMeta table. 
% INPUT
% key  = A key or other restriction of ns.Experiments to process only those. 
%        By default only those that don't alreay have info/nrframes metadata
%        will be affected. [struct]
% PV pairs:
% update = Force an update of experiment meta data.
arguments 
    key (1,1) = struct %  A restriction of experiments to consider. By default all experiments are processed.
    pv.update (1,1) logical = false  % Set to true to update the meta data
end

METANAME = {'nrframes';'xscale';'yscale'}; % These meta will be populated  (as string)
% Look for .sbx file in the ns.File table
sbxFiles = ((ns.File & key) & 'extension=''.sbx''');
if ~isstruct(key) % allow key to be a dj.Relvar like ns.Experiment
    key = fetch(key);
end

expts = (ns.Experiment & key ) & sbxFiles; 
if pv.update
    % Remove the info/nrframes meta data from the experiments
    % Then the code below will replace it.
    meta = (ns.ExperimentMeta & expts) & struct('meta_name',METANAME);
    delQuick(meta)
end

% Consider only files that don't have the meta data already
expts = expts - (ns.ExperimentMeta & struct('meta_name',METANAME));
for e = fetch(expts)'
   info = sbx.readInfoFile(e);   
   N= sbx.nrFrames(e,info);
    switch info.scanbox_version
        case 3
            xscale = info.dxcal;
            yscale = info.dycal;
            % % Calibration contains curvefit objects, which mym does not like. 
            % info.calibration.formulaX = formula(info.calibration.fx);
            % info.calibration.formulaY = formula(info.calibration.fy);
            % info.calibration.fx = [info.calibration.fx.a info.calibration.fx.b info.calibration.fx.c info.calibration.fx.d];            
            % info.calibration.fy = [info.calibration.fy.a info.calibration.fy.b info.calibration.fy.c info.calibration.fy.d];
            % 
        otherwise
            cal = info.calibration(info.config.magnification);
            xscale = 1./cal.x; 
            yscale = 1./cal.y;
    end
    % Insert in the ExperimentMeta table (as string)
    tpl = struct('subject',e.subject, ...
                    'session_date',e.session_date,...
                    'starttime',e.starttime, ...
                    'meta_name',METANAME, ...
                    'meta_value',{string(N),string(xscale),string(yscale)}');
    insert(ns.ExperimentMeta,tpl);
end
end