function addExperimentMeta(key,pv)
% After running nsScan to populate the ns.Experiments and ns.File table,
% this funciton will add the metadata from the scanbox .mat file (i.e. the 
% info struct plus the number of frames (nrframes) in the sbx file) to the
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

METANAME = {'info';'nrframes';'xscale';'yscale'}; % These meta will be populated
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

% Consider only files that don't have the info meta data 
expts = expts - (ns.ExperimentMeta & struct('meta_name',METANAME));
for e = fetch(expts)'
    fname = fetch1(sbxFiles & e,'filename');
    fldr = folder(ns.Experiment & e);
    ff =fldr + fname(1:end-3) + 'mat';
    if exist(ff,'file')
        load(ff ,'info'); % Get the info struct
    else
        error('File %s not found',strrep(ff,'\','/'))
    end
    N= sbx.nrFrames(fldr + fname,info);

    cal = info.calibration(info.config.magnification);
    xscale = 1./cal.x; 
    yscale = 1./cal.y;
    % Insert in the table
    tpl = struct('subject',e.subject, ...
                    'session_date',e.session_date,...
                    'starttime',e.starttime, ...
                    'meta_name',METANAME, ...
                    'meta_value',{info,N,xscale,yscale}');
    insert(ns.ExperimentMeta,tpl);
end
end