%{
# Subject information
subject: varchar(50)            # Code or number that identifies a subject 
---
sex = NULL : varchar(100)          # biological sex
dob = NULL : date                  # date of birth (ISO 8601)
species = NULL : varchar(100)      # Species
%}
% BK  - April 2022
classdef Subject < dj.Manual

    methods (Access = ?ns.Experiment)
        function subject = nwb(tbl,nwbRoot,pv)      
            % nwb export function, called from ns.Experiment nwbExport.      
            arguments
                tbl (1,1) ns.Subject {mustHaveRows(tbl,1)} %  A single row in this table
                nwbRoot (1,1) NwbFile  % The root element 
                pv (1,1) struct % The struct of pv pairs set in ns.Experiment/nwbExport
            end
            % Return the table as a NWB type. 
            tpl =fetch(tbl,'*');
           
            sexMapper = dictionary("MALE",'M',"FEMALE",'F',"UNKNOWN",'U', "OTHER",'O');
            if isKey(sexMapper,upper(tpl.sex))
                sex = sexMapper(upper(tpl.sex));
            else
                sex =upper(tpl.sex);
            end
            speciesMapper = dictionary("MOUSE",'Mus musculus',"HUMAN",'Homo sapiens');
            if isKey(speciesMapper,upper(tpl.species))
                species = speciesMapper(upper(tpl.species));
            else
                species =tpl.species;
            end
            
            subject = types.core.Subject( ...
                'subject_id', tpl.subject, ...
                'date_of_birth',datetime(tpl.dob,"TimeZone",pv.tz,"Format","uuuu-MM-dd ZZZZ","InputFormat","uuuu-MM-dd"),...
                'species', char(species), ...
                'sex', char(sex), ...
                'description',tpl.subject);
            
            %Loop over the subject meta dictionary
            for m = pv.subjectMeta.keys
                thisMeta =ns.getMeta(tbl,m);
                subject.(pv.subjectMeta(m)) =char(thisMeta.(m));
            end
            nwbRoot.general_subject = subject;
        end
    end
end