%{
# A table with the names of the data files for each experiment
-> ns.Experiment         # The experiment to which this belongs (FK)
filename: varchar(255)   # The relative filename
---
extension : varchar(10)  # File extension for easy filtering.
%}
%
% BK = April 2022

classdef File < dj.Computed

    methods (Access= protected)

        function makeTuples(tbl,key)
            % TODO handle p.Results.folderFun
            
            root = get(ns.Global,'root');
            % Search for files with matching prefix (subject.paradigm.startTime.*) in
            % the same folder and all files in the folder with this name.
            exptTpl = fetch(ns.Experiment &key,'*');
            pth = fullfile(root,strrep(exptTpl.session_date,'-','/'));
            prefix  = fullfile(pth,regexprep(exptTpl.file,'(\.mat$)','*')); % Swap extension
            inFolder = dir(prefix);
            % Search for a folder with matching prefix; add its content
            prefix = fullfile(prefix,'*');
            inSubFolder = dir(prefix);
            linkedFiles = cat(1,inFolder,inSubFolder);

            % Remove folders
            linkedFiles([linkedFiles.isdir]) =[];
            % Extract extension 
            match = regexp({linkedFiles.name},'.*\.(?<ext>.*$)','names');
            match = [match{:}];
            ext = strcat('.',{match.ext});
            
            % Add each one
            for f=1:numel(linkedFiles)                
                % Remove the part of the path that points to the folder
                % with the neurostim file, but keep subfolders deeper than
                % that (for most files relFolder will be '').
                relFolder = strrep(linkedFiles(f).folder,pth,'');
                qry = mergestruct(key,struct('filename',fullfile(relFolder,linkedFiles(f).name)));
                thisFile = ns.File & qry;
                if ~thisFile.exists
                    qry.extension = ext{f};
                    insert(tbl,qry);
                end
            end
        end
    end

end