%{
# A segmentation of the scan data obtained in a single session
-> ns.Session
-> sbx.Preprocessing
---
folder : varchar(1024) # Folder with the results of the preprocessing
img    : longblob       # Mean image
nrframes : int          # Total number of frames in the session.
%}

classdef Segmentation < dj.Computed

    properties (Dependent)
        keySource
    end
    methods
        function v = get.keySource(~)
            % Restrict to sessions that have sbx file
            v =(ns.Session & sbx.Scan)*sbx.Preprocessing;
        end
    end
    methods (Access=public)
        function v = getFolder(tbl)
             % subject.suite2p.name-of-preprocessing
             sessionPath=unique(folder(ns.Experiment & tbl));
             prep = fetch(sbx.Preprocessing & tbl,'name','toolbox');                
             session =fetch(ns.Session & tbl,'subject');
             v = fullfile(sessionPath, [session.subject '.' prep.toolbox '.' prep.name]);
             if ~exist(v,'dir')
                 fprintf('Folder %s  does not exist. Is NS_ROOT set correctly?\n',v);
             end
        end
    end
    methods (Access=protected)
        
        function makeTuples(tbl,key)                
                sessionPath=unique(folder(ns.Experiment & key));
                prep = fetch(sbx.Preprocessing & key,'*');                
                 % Set the output folder to be
                 % subject.suite2p.preprocessing 
                 % in the session folder
                resultsFolder = [key.subject '.' prep.toolbox '.' prep.name];
                % Find Experiments in this session that have Scans and
                % extract the folder name (subfolder named after the
                % Experiment).
                dataFldr = file(ns.Experiment & (sbx.Scan & key));                
                dataFldr = cellstr(strrep(dataFldr,'.mat',filesep))'; % cellstr to make py.list
                % Check that all folders exist.
                noDir = cellfun(@(x)exist(x,'dir'),dataFldr)==0;
                if any(noDir)
                    dataFldr{noDir} %#ok<NOPRT> 
                    error('SBX file folder not found');
                end
                switch (prep.toolbox)
                    case 'suite2p'
                        ops = py.suite2p.default_ops();
                        ops{'input_format'} = "sbx";
                        if prep.name ~= "default"                            
                            %replace parameters defined in the prep
                            %settings
                            % TODO
                        end
                        resultsFile =fullfile(sessionPath,resultsFolder,'plane0','ops.npy');
                        if ~exist(resultsFile,'file')
                        % Create a dict with the relevant information
                        db= py.dict(pyargs('save_path0',sessionPath, ...
                                            'save_folder',resultsFolder, ...
                                            'data_path',py.list(dataFldr), ...
                                            'fast_disk',fullfile(sessionPath,resultsFolder)));
                        % Pass to python to process
                        fprintf('Starting suite2p run_s2p at %s... this will take a while \n',datetime('now'))
                        py.suite2p.run_s2p(ops =ops,db=db);
                        fprintf('Completed at %s\n',datetime('now'));
                        else
                            fprintf('Segmentation results already exists. Importing %s\n',resultsFile);
                        end
                        % Load the save ops.npy to extract the mean image
                        ops =py.numpy.load(resultsFile,allow_pickle=true);
                        img= single(ops.item{'meanImg'}); % Convert to single to store in DJ
                        N = double(ops.item{'nframes'});
                        tpl = mergestruct(key,struct('img',img,'folder',fullfile(sessionPath,resultsFolder),'nrframes',N));
                        insert(tbl,tpl);
                    case 'caiman'
                        % TODO
                    otherwise
                        error('Unknown preprocessing toolbox %s',prep.toolbox);
                end


        end
    end
end