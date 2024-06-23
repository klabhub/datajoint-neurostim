function nwbTrg = createNwb(expt,pv)
% Create Neurodata without Borders objects for each experiment in a
% ns.Experiment table. 
% 
% folder = Folder where the data files will be saved. By default this is
% empty, which will create and return the NwbFile objects, but not save
% them.
% tz  - TimeZone to use for session_Start_time and
% timestamps_reference_time ["local"]
% general -  A struct with fields named after the general_ properties of the
% NwbFile object that use free form text to describe meta data. 
% The allowable names depend on the NWM file schema.  By default this is empty. 
%
%
% EXAMPLE

%  general = struct('lab',"KLab",'Insitution',"Rutgers University - Newark")
% 
% The columns of the ns.Subject table are always included, to include ns.SubjectMeta
% data as well, provide a dictionary that maps the names of the meta data
% to names in the NWB subject schema.
% 
% subjectMeta =dictionary('strain','strain;')
% 
% To include raw data , provide a dictionary from file extensions (as used
% in ns.File) to function handles that take and NwbFile object and the
% full filename with the raw data as their input.  See sbx.nwb for an
% example of such a function. Files in ns.File that do not match an
% xtension in the raw dictionary will be ignored. If an experiment does not
% have a corresponding raw data file in ns.File, a warning will be issued.
% 
%  raw = dictionary("sbx",@sbx.nwb) ;
%
% createNwb(ns.Experiment ,general=general, subjectMeta= subjectMeta,
%                           raw=raw, folder = "c:/temp");

arguments  
    expt (1,1) ns.Experiment
    pv.folder (1,1) = "" 
    pv.tz (1,1) string = "local"
    pv.general (1,1) struct = struct();    
    pv.subjectMeta (1,1) dictionary  = dictionary(string([]),string([]));
    pv.raw (1,:) dictionary =dictionary(string([]),@(x) (x)) 
end

assert(~isempty(which('NwbFile')),'This function depends on the matnwb package. Install it from github and add it to the Matlab path');


for e = fetch(expt,'*')'
    % Loop over experiments (NWB refers to this as a "session", NS uses session to refer to all experiments for 
    % a subject on a given day)
    uniqueExperimentName = sprintf('%s_%s_%s',e.subject,e.session_date,e.starttime);
    
    
    nwbTrg = NwbFile(...
        'session_description',sprintf('%s',e.paradigm ),...
        'identifier',char(java.util.UUID.randomUUID().toString()),...
        'session_start_time', datetime(e.session_date,'TimeZone',pv.tz),...
        'timestamps_reference_time',datetime([e.session_date 'T' e.starttime],'TimeZone',pv.tz), ...        
        'general_session_id', uniqueExperimentName);

    
    
        %% General properties (specified in the call to this function)
        fn = fieldnames(pv.general);
        for i=1:numel(fn)
            prop = "general_" + fn{i};
            try
                nwbTrg.(prop) = pv.general.(fn{i});
            catch  
                fprintf(2,"Property %s does not exist in the NWB schema. Ignored. \n", prop );
            end
        end

        %% Raw data,read from ns.File and matched against the user provided dictionary.
        fldr = folder(expt&e);
        for rawExt = keys(pv.raw)
            f =fetch(ns.File & e & struct('extension',rawExt));
            if isempty(f)
                fprintf("No %s raw data found in %s \n" ,rawExt,uniqueExperimentName);
            else
                % Call the user provided function to update the nwbFile
                % object
                fun = pv.raw(rawExt);
                nwbTrg =fun(nwbTrg,fullfile(fldr,f.filename));
            end
        end


        if ~isempty(pv.folder)
            assert(exist(pv.folder,"dir"),"Folder %s does not exist.",pv.folder);
            fname = fullfile(pv.folder,[strrep(uniqueExperimentName,':','') '.nwb']);
            nwbExport(nwbTrg,fname);
        end

end
