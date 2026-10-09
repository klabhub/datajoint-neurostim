%{
# Fiducial locations in voxel coordinates
-> ns.Subject   
mri_date : date  #Date of the MRI acquisition (ISO 8601)
coordsys : varchar(10) # E.g. ctc for nas,lpa,rpa.
---
fiducial: blob  # Struct with voxel coordinates of the fiducials
folder : varchar(255) # Session folder that stores the dicoms
%}

classdef Fiducial <dj.Manual

    methods (Static)
        function define(pv)
            arguments
                pv.qry (1,1) = [] % DJ query that returns subjects - overrules pv.subject
                pv.subject (1,:) = string.empty
                pv.newOnly (1,1) logical = true
                pv.dicomSubFolder (1,:) string = ["MPRage" "t1_mpr"] % Used to find subfolders 
                pv.coordsys (1,1) string = "ctf"
            end

            assert(exist("ft_read_mri","file"),"Please add FieldTrip to your path.")

            if isa(pv.qry,"dj.Relvar")
                % Use a query to define the subjects
                assert(isempty(pv.subject),"subject should be empty when defining a qry");
                pv.subject = string(fetchn(pv.qry,'subject'))';
            end

            %% Find candidate folders with MRI data
            % (they are named subject.dicoms)
            if isempty(pv.subject)
                % Select all
                subjects =fetchtable(ns.Subject,'subject');
                pv.subject = unique(subjects.subject)';
            end
            if pv.newOnly
                subjectsWithFiducials =fetchtable(sloc.Fiducial,'subject');
                pv.subject= setdiff(pv.subject,subjectsWithFiducials.subject);
            end
            if isempty(pv.subject)
                fprintf('No subjects remaining after newOnly\');
                return;
            end
            
            % Loop over remaining subjects
            for subject = pv.subject
                sessions = ns.Session & struct('subject',subject);
                fldrs = fullfile(folder(sessions) ,subject+ ".dicoms");
                keys = fetch(sessions);
                nrFldrs = numel(fldrs);
                stay = false(nrFldrs,1);
                for d=1:nrFldrs
                    stay(d) = exist(fldrs(d),"dir");
                end
                fldrs(~stay) = [];
                keys(~stay) = [];

                nrFldrs = numel(fldrs);
                if isempty(fldrs)
                    warning('No DICOM folders found.');
                    continue;
                end

                for i = 1:nrFldrs
                    % read the mri
                    fprintf('Working on %s \n',fldrs(i));
                    % save it to .nii for future faster reads
                    niftiFile = fullfile(fileparts(fldrs(i)),subject + "_anat.nii");
                    dcmFile  = sloc.Fiducial.firstDicomFile(fldrs(i),pv.dicomSubFolder);
                    if isempty(dcmFile)
                        fprintf(2,'No dicom files found in %s. Dicom subfolder mismatch? Skipping.\n',fldrs(i));
                        continue;
                    end
                    if exist(niftiFile,"file")
                        fprintf('Reading from nifti file %s\n',niftiFile)
                        mri = ft_read_mri(char(niftiFile));
                    else
                        mri = sloc.Fiducial.dicom2nifti(fullfile(dcmFile.folder,dcmFile.name),niftiFile);                        
                    end
                    %Always get the study date from the first dicom file
                    %(not in the nifti?).
                    info = dicominfo(fullfile(dcmFile.folder,dcmFile.name));
                    studyDate = datetime(info.StudyDate,InputFormat= 'uuuuMMdd',Format='uuuu-MM-dd');


                    % Use FT interactive mode to determine fiducials
                    cfg = [];
                    cfg.method = 'interactive';
                    cfg.coordsys = char(pv.coordsys);% the desired coordinate system
                    mri_realigned = ft_volumerealign(cfg, mri);


                    % Construct tuple to add to the table.
                    fiducial = mri_realigned.cfg.fiducial;
                    for fn = string(fieldnames(fiducial))'
                        if all(isnan(fiducial.(fn)))
                            fiducial = rmfield(fiducial,fn);
                        end
                    end
                    if isempty(fieldnames(fiducial))
                        fprintf('No fiducials left. Skipping %s\n ',niftiFile);
                    else
                        thisKey  = keys(i);
                        thisKey  = rmfield(thisKey,'session_date');
                        thisKey.folder = strrep(extractAfter(fldrs(i),getenv('NS_ROOT')),'\','/');
                        thisKey.mri_date = char(studyDate);
                        thisKey.fiducial = fiducial;
                        thisKey.coordsys = cfg.coordsys;
                        insert(sloc.Fiducial,thisKey);
                    end
                end
            end
        end
    end
    methods (Static, Access=public)
        function mri = dicom2nifti(firstDicomFile,niftiFile)
            % Read the dicoms and convert to nifti
            mri = ft_read_mri(firstDicomFile, 'dataformat','dicom');
            mri = ft_convert_units(mri,'mm');
            % Save to nifti for future faster reads 
            ft_write_mri(char(niftiFile),mri,'dataformat', 'nifti');                        
        end

        function first = firstDicomFile(fldr,subFolders)
            arguments
                fldr (1,1) string
                subFolders (1,:) string = ["MPRage" "t1_mpr"] % Used to find subfolders with contains
            end
            % Find the dicom subfolder in fldr that matches one of the subFolders
            % (case-insensitive)            
            candidates = dir(fldr);
            candidates = candidates([candidates.isdir]);
            isCandidate = contains({candidates.name},subFolders,'IgnoreCase',true);
            candidates = candidates(isCandidate);
            assert(isscalar(candidates),"No unique match for a dicom subfolder in %s",fldr);
            dcmFiles = dir(fullfile(fldr,candidates.name,"*.dcm"));    
            if isempty(dcmFiles)
                fprintf(2,'No dicom files found in %s. Dicom subfolder mismatch? Skipping.\n',fldr);
                first = struct([]);
            else
                first = dcmFiles(1);               
            end
        end
    end
end