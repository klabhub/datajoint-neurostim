function str = limitDate(dt,schedule)
% Given a date and a schedule, this function returns the MySql statement
% that selects the matching sessions.
%
% dt - a datetime value
% schedule - "d","m","y","w" to select the day, month, year, or week that
%               the dt is part of.
%           For "w" the week is defined as the preceding 7 days.
%           "a"  - All times 
% BK - Oct 2023.
arguments
    dt (1,1) datetime
    schedule (1,1) string {mustBeMember(schedule,["d" "m" "y" "w"])}
end

fmt = "yyyy-MM-dd";
dt.Format = fmt;

switch (schedule)
    case "a"
        % All 
        startDate = datetime("1900-01-01");
        stopDate = tdatetime("tomorrow");
    case "d"
        startDate = dt;
        stopDate = dt;
    case "w"
        startDate = dt-7
        stopDate =dt;
    case "m"
        startDate = datetime(year(dt),month(dt),1);
        stopDate =  datetime(year(dt),month(dt),31);
    case "y"
        startDate  =datetime(year(dt),1,1);
        stopDate  = datetime(year(dt),12,31);
end

startDate.Format  =fmt;
stopDate.Format  =fmt;
str = ['session_date BETWEEN ''' char(startDate) ''' AND  ''' char(stopDate) ''''];
