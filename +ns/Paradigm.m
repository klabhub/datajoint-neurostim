%{
#  Paradigms associated with this database
name         :  varchar(32)     # A  unique name for this paradigm
---
description = NULL : varchar (255)     # Short description
mintrials  = 1  : smallint          # Only experiments with at least this number of trials are valid
from = NULL : datetime          # Only experiments from this day onward (optional)
to   = NULL : datetime          # Only experiments until this day (optional)
%}
% This table is used by nsScan to determine whether a Neurostim experiment
% should be included or not. Required selection is based on the paradigm
% name but additional inclusion criteria can be added by specifying
% mintrials (the minimum number of trials in the experiment), from (earliest
% date), to (last date).
% For instance, to include experiments with paradigm 'xxx' and at least 10
% trials (but on any date).
% insert(ns.Paradigm,struct('name','xxx','mintrials',10))
% Note that paradigm names are matched case-insensitively and that the paradigm name is 
% truncated (with a warning) to 32 characters if longer.

classdef Paradigm < dj.Lookup
    
    methods (Access = public )
        function insert(tbl, tpl)
            % Check the length of the name; to avoid problems with PK length for tables that link to this or the
            % experiment table, paradigms are 32 characters long. If the name is longer than 32 characters, we truncate it and issue a warning.
           if length(tpl.name) > 32
                    warning('Paradigm name is longer than 32 characters. It will be truncated.');
           end
           tpl.name = tpl.name(1:min(32,length(tpl.name)));
           insertIfNew(tbl,makeMymSafe(tpl));           
        end
    end


    methods (Static)
        function tbl = mismatch()
            % Find ns.Experiment entries that should not have been added if
            % the current ns.Paradigm table restrictions had been applied
            % to nsScan always
            qry ='';
            for f= fetch(ns.Paradigm,'*')'
                if ~isempty(f.from)
                    frmQry = sprintf(' | session_date < %s',f.from) ;
                else 
                    frmQry = '';
                end
                                    
                if ~isempty(f.to)
                    toQry = sprintf('| session_date > %s',f.to) ;
                else
                    toQry ='';
                end
                thisQry = sprintf('(paradigm="%s" & (nrtrials < %d %s %s))',f.name,f.mintrials,frmQry,toQry);
                qry = [qry  thisQry '|'];
            end
            qry(end)=[];
            tbl = ns.Experiment & qry;
        end
    end

end
