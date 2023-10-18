%{
# A table with the names of the data files for each experiment
-> ns.Experiment         # The experiment to which this belongs (FK)
filename: varchar(255)   # The relative filename (using / to separate folders)
---
extension : varchar(10)  # File extension for easy filtering.
%}
%
% BK = April 2022

classdef File < dj.Imported

    methods (Access = public)
        
        function out = checkExists(tbl)
            % Check whether files exist on the current file system
            root = getenv('NS_ROOT');
            if isempty(root)
                fprintf('NS_ROOT is empty. Most likely no files will be found.\n')
            end            
            nrFiles = count(tbl);
            exists = false(nrFiles,1);
            filename =repmat("",[nrFiles 1]);
            T=table(filename ,exists);
            fCntr=0;
            for f = tbl.fetch()'
                fCntr =fCntr+1;
                fldr = folder(ns.Experiment &f);
                full = fullfile(fldr,f.filename);
                T.filename(fCntr)= full;
                T.exists(fCntr)= exist(full,"file")==2;
            end
            if nargout==0
                groupcounts(T,"exists")
            else
                out =T;
            end
        end
    end
    methods (Access= protected)

        function makeTuples(tbl,key)
            % TODO handle p.Results.folderFun

            root = getenv('NS_ROOT');
            % Search for files with matching prefix (subject.paradigm.startTime.*) in
            % the same folder and all files in the folder with this name.
            exptTpl = fetch(ns.Experiment &key,'*');
            pth = fullfile(root,strrep(exptTpl.session_date,'-','/'));
            prefix  = fullfile(pth,regexprep(exptTpl.file,'(\.mat$)','*')); % Swap extension
            inFolder = dir(prefix);
            % Search for a folder with matching prefix; add its content
            % (inluding content in all sub/sub folders).
            subFolder = inFolder([inFolder.isdir]);
            inSubFolder = dir(fullfile(subFolder.folder,subFolder.name,'**','*'));
            linkedFiles = cat(1,inFolder,inSubFolder);
            % Remove folders
            linkedFiles([linkedFiles.isdir]) =[];
            if ~isempty(linkedFiles)
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
                    filename  = strrep(fullfile(relFolder,linkedFiles(f).name),'\','/'); % Force / convention.
                    qry = mergestruct(key,struct('filename',filename));
                    thisFile = ns.File & qry;
                    if ~thisFile.exists
                        qry.extension = ext{f};
                        insert(tbl,qry);
                    end
                end
            end
        end
    end

end