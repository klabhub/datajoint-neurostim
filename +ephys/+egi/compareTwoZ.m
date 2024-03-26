function [zBeforeFile, zAfterFile] = compareTwoZ(mffInfo,s,stay,zBeforeFile,mffFile)
% Since we only have one zCheck, the mffFile likely contains an impedance
% check. Find out if the zCheck was performed before or after the
% experiment and correct 'zAfterFile' and 'zBeforeFile' accordingly
if nargin <5
    sibwarn('File with EGI plugin, but no matching MFF file found. Candidate MFF files do not match in terms of day/paradigm');
    char(mffInfo.name)
end

    neurostimTime =  datenum([s.date ' ' s.time] ,'dd-mmm-yyyy HH:MM:SS');
    
    try
        tempTimeString = mffInfo(contains({mffInfo.name},strrep(zBeforeFile,[fileparts(zBeforeFile) '\'],''))).time; %local
    catch
     	tempTimeString = mffInfo(contains({mffInfo.name},strrep(zBeforeFile,[fileparts(zBeforeFile) '/'],''))).time; %cluster
    end
    
    zCheckTime = [tempTimeString(1:2) ':' tempTimeString(3:4) ':' tempTimeString(5:6)];
    zCheckTime = datenum(zCheckTime,'HH:MM:SS');
    
if zCheckTime<neurostimTime
    zAfterFile = mffFile;
else
    zAfterFile = zBeforeFile;
    zBeforeFile = mffFile;
end

end
