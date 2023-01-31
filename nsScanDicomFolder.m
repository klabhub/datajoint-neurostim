function meta = nsScanDicomFolder(folderFullName)
% Example function that receives a list of folders and returns a struct
% array with information on the experiments in that folder. This is used,
% for instance, to include experiments that did not have any real neurostim
% experiment associated with them (e.g., a structural MRI).
%
%
% To be used by nsScan, the struct array must have the following fields:
% .subject
% .starttime
% .session_date
% .paradigm  
% .isFolder = true  
% .type    (some indicator of the file type)
% .bytes  (size of the file)
% .provenance (history of file name changes).
% Dec 2022
folderExtension = '.dicoms'; % This function will search for folders with XXX.dicoms format where xxx is the subject iD
meta =struct([]);
if strcmpi(filesep','\')
    fs = '\\';
else
    fs = filesep;
end
folderPattern = ['(?<root>.*)\\(?<session_date>\d{4,4}' fs '\d{2,2}' fs '\d{2,2})' fs '(?<subject>\w{2,10})\'  folderExtension '$'];
dataFolders = regexp(folderFullName,folderPattern,'names');
out = cellfun(@isempty,dataFolders);
dataFolders = [dataFolders{~out}];
for i=1:numel(dataFolders)
    subs = dir(fullfile(dataFolders(i).root,dataFolders(i).session_date,[dataFolders(i).subject folderExtension],'*'));
    subs(~[subs.isdir] | ismember({subs.name},{'.','..'}))=[];
    for j=1:numel(subs)
        thisStruct = dataFolders(i);
        allFiles= dir(fullfile(subs(j).folder,subs(j).name,'*.dcm'));
        if ~isempty(allFiles)
            % Read the first one to extract start time
            info=dicominfo(fullfile(subs(j).folder,subs(j).name,allFiles(1).name));
            % SeriesTime is used to identify startTime of the experiment. It is
            % possible that two scans start within 0.5 s of each other (e.g. a
            % field map), then this will generate identical starttimes, which
            % will be detected as non unique PK below. When that happens, need
            % to think of a way to make this unique
            thisStruct.starttime = num2str(round(str2double(info.SeriesTime)));
            thisStruct.paradigm = subs(j).name;
            thisStruct.isFolder = true;
            thisStruct.type = '.dicoms';
            thisStruct.bytes= subs(j).bytes;
            meta= cat(2,meta,thisStruct);
        end
    end
end



if ~isempty(meta)
    pk = strcat({meta.session_date}','-', {meta.subject}','-',{meta.starttime}');
    [upk,stayIx] = unique(pk);
    delta = numel(pk)- numel(upk);
    if delta ~=0
        fprintf(2, '%d non unique experiments  (probably a start time that differed by less than 0.5 s. The duplicates have been removed. This  requires manual intervention.',delta)
        meta = meta(stayIx);
    end
    
    % Not needed in caller.
    meta = rmfield(meta,'root');
end