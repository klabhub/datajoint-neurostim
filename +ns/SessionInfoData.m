% -------------------------------------------------------------------------
% ns.SessionInfoData
% -------------------------------------------------------------------------
%{
# ns.SessionInfoData: Session-level metadata
-> ns.SessionInfo
-> ns.C
---
session_no : int
session_day: int # order of the day of the session
time_since_prev: varchar(32) # duration since the prev session in terms of 'HH:MM:SS' 
%}
classdef SessionInfoData < dj.Part & dj.DJInstance
    properties(SetAccess=protected)
        master = ns.SessionInfo
    end
end
