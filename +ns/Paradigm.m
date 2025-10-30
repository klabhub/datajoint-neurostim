%{
#  Paradigms associated with this database
name         :  varchar(255)     # A  unique name for this paradigm
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
% Note that paradigm names are matched case-insensitively.


classdef Paradigm < dj.Lookup

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
