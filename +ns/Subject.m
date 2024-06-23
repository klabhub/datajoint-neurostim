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

    methods (Access = public)

        function subject = nwb(tbl,meta)      
            arguments
                tbl (1,1) ns.Subject
                meta (1,1) dictionary = dictionary(string([]),string([]))
            end
            % Return the table as a NWB type. 
            sCntr= 0;
            for s =fetch(tbl,'*')'
                sCntr= sCntr+1;
            subject(sCntr) = types.core.Subject( ...
                'subject_id', s.subject, ...
                'date_of_birth',s.dob,...
                'age_reference','birth',...                                
                'species', s.species, ...
                'sex', s.sex); %#ok<AGROW>
            end
            for m = meta.keys   

            end
        end
    end
end