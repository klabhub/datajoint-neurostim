function addMeta(trg,info,pv)
% Function to add sbx specific meta data to the experiments table.
% Called automatically from sbx.Preprocessed/makeTuples when data
% are preprocessed (but can be called manually by calling this function
% with a Session, Experiment, or sbx.Preprocessed table. Note that the meta
% data are always added for all analyzeable experiments in the sessions in
% the input table.
% By default only new meta data are added, experiments with existing
% nrframes metadata are skipped.
%
% 
arguments
    trg  % A session, experiment, or sbx.Preprocessed table, or a single experiment tpl.
    info struct {mustBeScalarOrEmpty} = struct([]);
    pv.newOnly (1,1) logical = true  % Only add new tuples
    pv.verbose (1,1) logical = false 
end


if isempty(info) && isa(trg,'dj.Relvar')    
    trg  = ns.Session & trg; % Get all sessions
    % Limit to experiments with sbx files 
    experimentsWithSbx  = ns.Experiment & (ns.File & 'extension=".sbx"');
    trg = trg & (experimentsWithSbx);
    if pv.newOnly
        % Restrict to experiments that do not have the meta data yet.
        experimentsWithout = experimentsWithSbx - (ns.ExperimentMeta & 'meta_name="nrframes"');
        trg = trg & experimentsWithout;
    end
        
    for s = fetch(trg)' % Loop over sessions that have an entry in the table
        allExptThisSession = ns.Experiment & (ns.File & 'extension=''.sbx''') &s;
        analyzeExptThisSession = analyze(allExptThisSession,strict=false);
        for e=fetch(analyzeExptThisSession)' % Loop over experiments in the session
            info        = sbx.readInfoFile(e); % Read info file
            if pv.verbose
                fprintf('Adding SBX meta data for %s on %s at %s\n',e.subject,e.session_date,e.starttime);
            end
            addMetaForExpt(e,info);
        end
    end
elseif isstruct(trg)
    if pv.newOnly && exists( (ns.ExperimentMeta & 'meta_name="nrframes"' & trg))
        return;
    else
        addMetaForExpt(trg,info);
    end
else
    error('sbx.addMeta input arguments not correct')
end
end

function  addMetaForExpt(expt,info)
% Store the nrframes and nrplanes and depth as meta data for quicker
% access and sanity checks.
arguments
    expt (1,1) struct
    info (1,1) struct
end

metaE = ns.ExperimentMeta & expt & struct('meta_name',{'nrframes','nrplanes','depth'}');
tpl = mergestruct(expt,struct('meta_name',{'nrframes','nrplanes','depth'}','meta_value',{num2str(info.nrFrames),num2str(info.nrPlanes),num2str(info.config.knobby.pos.z)}'));

dj.conn().startTransaction
if exists(metaE)
    delQuick(metaE); % Replace below
end
insert(ns.ExperimentMeta,tpl);
dj.conn().commitTransaction;
% Because meta data are assumed to originate from json (and are
% sometimes deleted to be replaced with values from json)
% save to json too so that any del/replace operation will work.
jsonFile = strrep(file(ns.Experiment &expt),'.mat','.json');
if exist(jsonFile,"file")
    metaData = readJson(jsonFile);
    metaData(1).nrframes= info.nrFrames; % use (1) in case the metaDat read from file is empty struct
    metaData(1).nrplanes = info.nrPlanes;
    metaData(1).depth =info.config.knobby.pos.z;
else
    metaData = struct('nrframes',info.nrFrames,'nrplanes',info.nrPlanes,'depth',info.config.knobby.pos.z);
end
writeJson(metaData,jsonFile);
fprintf("Wrote %s\n",jsonFile);
end

