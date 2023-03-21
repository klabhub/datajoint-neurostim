%{
# A complete preprocessed data set of the SBX data obtained in a single session.
-> ns.Session
-> sbx.PrepParms
---
folder : varchar(1024)      # Folder with the preprocessing results
img    : longblob           # Mean image
nrframesinsession : int     # Total number of frames in the session.
%}
%
% Sessions with entries in the Populate this table
classdef Preprocessed < dj.Imported
    properties (Dependent)
        keySource
        ops
        stat
    end
    methods
        function v = get.keySource(~)
            % Restrict to sessions that have sbx file
            v =(ns.Session & (ns.File & 'extension=''.sbx'''))*sbx.PrepParms;
        end

        function v =  get.ops(tbl)
            keyCntr = 0;            
            for key =fetch(tbl*sbx.PrepParms,'*')
                  keyCntr = keyCntr+1;
             sessionPath=unique(folder(ns.Experiment & key));
             resultsFolder = [key.subject '.' key.toolbox '.' key.prep];
             resultsFile =fullfile(sessionPath,resultsFolder,'plane0','ops.npy');            
              v{keyCntr} =   py.numpy.load(resultsFile,allow_pickle=true); %#ok<AGROW>              
            end
            if numel(v)==1
                v =v{1};
            end
                 
        end

        function v =  get.stat(tbl)
            keyCntr = 0;
            for key =fetch(tbl*sbx.PrepParms,'*')
              keyCntr = keyCntr+1;
              sessionPath=unique(folder(ns.Experiment & key));
             resultsFolder = [key.subject '.' key.toolbox '.' key.prep];
             resultsFile =fullfile(sessionPath,resultsFolder,'plane0','stat.mat');
              load(resultsFile,'stat');
              v{keyCntr} =   [stat{:}];            %#ok<PROP> 
            end
            if numel(v)==1
                v =v{1};
            end
        end

    end
    methods (Access=public)
        function v = getFolder(tbl)
            % subject.suite2p.name-of-preprocessing
            sessionPath=unique(folder(ns.Experiment & tbl));
            prep = fetch(sbx.PrepParms & tbl,'prep','toolbox');
            session =fetch(ns.Session & tbl,'subject');
            v = fullfile(sessionPath, [session.subject '.' prep.toolbox '.' prep.prep]);
            if ~exist(v,'dir')
                fprintf('Folder %s  does not exist. Is NS_ROOT set correctly?\n',v);
            end
        end

        
    end
    methods (Access=protected)

        function makeTuples(tbl,key)
            sessionPath=unique(folder(ns.Experiment & key));
            prep = fetch(sbx.PrepParms & key,'*');
            % Set the output folder to be
            % subject.suite2p.preprocessing
            % in the session folder
            resultsFolder = [key.subject '.' prep.toolbox '.' prep.prep];
            % Find Experiments in this session that have Scans and
            % extract the folder name (subfolder named after the
            % Experiment).
            dataFldr = file(ns.Experiment & (ns.File & 'extension=''.sbx''' & key));
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
                        %replace parameters defined in the prep
                        %settings
                        fn= fieldnames(prep.parms);
                        for  f= 1:numel(fn)
                            try
                                current= ops{fn{f}};
                                new = feval(class(current),prep.parms.(fn{f}));
                            catch me
                                error('Parameter %s does not exist in default_ops(). Typo?',fn{f});
                            end
                            ops{fn{f}} = new;
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
                        % Couldn't figure out how to convert stat.npy so
                        % save it as .mat (Not using save_mat to avoid
                        % duplicating all of the fluorescence data in F.mat
                        % and I don't want to delete the .npy files because
                        % they are useful to view in the suite2p gui. 
                        statFile = fullfile(sessionPath,resultsFolder,'plane0','stat.npy');
                        npyToMat(statFile);
                        fprintf('Completed at %s\n',datetime('now'));
                    else
                        fprintf('Preprocessing results already exists. Importing %s\n',resultsFile);
                    end
                    % Load the save ops.npy to extract the mean image
                    ops =py.numpy.load(resultsFile,allow_pickle=true);
                    img= single(ops.item{'meanImg'}); % Convert to single to store in DJ
                    N = double(ops.item{'nframes'});
                    tpl = mergestruct(key,struct('img',img,'folder',fullfile(sessionPath,resultsFolder),'nrframesinsession',N));
                    insert(tbl,tpl);
                case 'caiman'
                    % TODO
                otherwise
                    error('Unknown preprocessing toolbox %s',prep.toolbox);
            end


        end
    end
    methods (Static)
    
    end
end