%{
# A table with the names of the data files for each experiment
-> ns.Experiment         # The experiment to which this belongs (FK)
filename: varchar(255)   # The relative filename (using / to separate folders)
---
extension : varchar(10)  # File extension for easy filtering.
bytes = NULL : float     # Bytes in the file
checksum = NULL : char(32) # MD5 Hash checksum
%}
%
% BK = April 2022

classdef File < dj.Imported

    methods (Access = public)

        function out = checkExists(tbl,pv)
            arguments
                tbl (1,1) ns.File
                pv.bytes (1,1) logical = false
                pv.checksum (1,1) logical = false
            end

            % Check whether files exist on the current file system
            root = getenv('NS_ROOT');
            if isempty(root)
                fprintf('NS_ROOT is empty. Most likely no files will be found.\n')
            end
            nrFiles = count(tbl);
            exists = false(nrFiles,1);
            filename =repmat("",[nrFiles 1]);
            T=table(filename ,exists);
            if pv.bytes
                T =addvars(T,false(nrFiles,1),'newVariableNames','bytes');
            end
            if pv.checksum
                T =addvars(T,false(nrFiles,1),'NewVariableNames','checksum');
            end
            fCntr=0;
            for f = tbl.fetch('bytes','checksum')'
                fCntr =fCntr+1;
                fldr = folder(ns.Experiment &f);
                full = fullfile(fldr,f.filename);
                T.filename(fCntr)= full;
                T.exists(fCntr)= exist(full,"file")==2;
                if T.exists(fCntr)
                    % File exists
                    if pv.bytes
                        % Check bytes ifrequested.
                        if isnan(f.bytes)

                        else
                            d = dir(full);
                            T.bytes(fCntr) = f.bytes==d.bytes;
                        end
                    end
                    if pv.checksum
                        % Check md5 checksum
                        if isnan(f.checksum)

                        else
                            md5 = ns.File.checksum(full);
                            T.checksum(fCntr) = string(f.checksum)==md5;
                        end
                    end
                end
            end
            if nargout==0
                groupcounts(T,"exists")
            else
                out =T;
            end
        end

        function updateWithFiles(tbl,key,linkedFiles)
            % Usually called from makeTuples to add the standard files
            % named according to NS conventions, but can be called to add
            % other file,s that don't match that format.
            linkedFiles([linkedFiles.isdir])=[];
            pth =  folder(ns.Experiment &key);
            for f=1:numel(linkedFiles)

                % Add byte counts and checksum
                ff = fullfile(linkedFiles(f).folder,linkedFiles(f).name);
                md5Hash = ns.File.checksum(ff);
                d = dir(ff);
                bytes = d.bytes;

                [~,~,ext] =fileparts(linkedFiles(f).name);
                % Remove the part of the path that points to the folder
                % with the neurostim file, but keep subfolders deeper than
                % that (for most files relFolder will be '').
                relFolder = strrep(linkedFiles(f).folder,pth,'');
                filename  = strrep(fullfile(relFolder,linkedFiles(f).name),'\','/'); % Force / convention.
                qry = mergestruct(key,struct('filename',filename));
                thisFile = ns.File & qry;
                if ~thisFile.exists
                    qry.extension = ext;
                    qry.bytes =bytes;
                    qry.checksum = md5Hash;
                    insert(tbl,qry);
                end
            end
        end



    end


    methods (Access= protected)

        function makeTuples(tbl,key)
            % TODO handle p.Results.folderFun

            % Search for files with matching prefix (subject.paradigm.startTime.*) in
            % the same folder and all files in the folder with this name.
            exptTpl = fetch(ns.Experiment &key,'*');
            pth =  folder(ns.Experiment &key);
            prefix  = fullfile(pth,regexprep(exptTpl.file,'(\.mat$)','*')); % Swap extension
            inFolder = dir(prefix);
            % Search for a folder with matching prefix; add its content
            % (inluding content in all sub/sub folders).
            subFolder = inFolder([inFolder.isdir]);
            linkedFiles= inFolder(~[inFolder.isdir]);
            for i=1:numel(subFolder)
                inSubFolder = dir(fullfile(subFolder(i).folder,subFolder(i).name,'**','*'));
                linkedFiles = cat(1,linkedFiles,inSubFolder(~[inSubFolder.isdir]));
            end
            updateWithFiles(tbl,key,linkedFiles);
        end
    end

    methods (Static)
        function md5Hash = checksum(ff)
            % Determine the MD5 checksum of a , using a JAVA library
            % Returns all zeros if the file does not exist.
            mdLibrary = java.security.MessageDigest.getInstance('MD5');
            if exist(ff,"file")
                d= dir(ff);
                if d.bytes>1e9
                    fprintf(2,"%s is bigger than 1 GB. Not computing MD5 checksum.\n",ff);
                    md5Hash = string(repmat('0',[1 32]));
                else
                    fid = fopen(ff,'r');
                    try
                        digest = dec2hex(typecast(mdLibrary.digest(fread(fid, inf, '*uint8')),'uint8'));
                        fclose(fid);
                        md5Hash = string(lower(digest(:)'));
                    catch me
                        fclose(fid);
                        fprintf(2,"MD5 failed on %s. (%s)\n",ff,me.message);
                        md5Hash = string(repmat('0',[1 32]));
                    end
                end
            else
                fprintf(2,"%s does not exist. Could not determine md5 hash\n",ff)
                md5Hash = string(repmat('0',[1 32]));
            end
        end
    end

end