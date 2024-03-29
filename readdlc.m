function [T,bodyparts,scorer] = readdlc(csvFile)
% function T = readdlc(csvFile)
% 
% Read the csv output file that contains pose information extracted by
% DeepLabCut and returns the contents as a table with columns named after
% the bodyparts. E.g.  if nose is tracked, T will have a column of nosex
% nosey and noselikelihood.
%
arguments
    csvFile (1,1) string 
end
if ~exist(csvFile,"file")
    error('%s not found',csvFile);
end

% Read the three header lines
fid =fopen(csvFile,'r');
scorer =strsplit(fgetl(fid),',');
scorer = unique(scorer(2:end));
if numel(scorer)==1
    scorer=scorer{1};
else
    scorer = strjoin(scorer,'/');
end
bodyparts = strsplit(fgetl(fid),',');
header = strsplit(fgetl(fid),',');
fclose(fid);
% Read the rest as a table
rawT = readmatrix(csvFile,'NumHeaderLines',3,'FileType','text','Delimiter',',');
isX = strcmpi(header,'x');
x = rawT(:,isX);
isY = strcmpi(header,'y');
y = rawT(:,isY);
isl = strcmpi(header,'likelihood');
likelihood = rawT(:,isl);
T =table(x,y,likelihood);
bodyparts =bodyparts(isX); % Just need them once
end