%{
# An experiment refers to a single run of a specific experimental paradigm
starttime: time # Time that the experiment started (HH:MM:SS)
-> ns.Session   # Corresponding session session  (FK)
---
paradigm: varchar(255)      # Name of the paradigm
file = NULL : varchar(255)  # File that contains the Neurostim output
stimuli = NULL : smallint      # Number of stimuli in the experiment
blocks = NULL : smallint       # Number of blocks in the experiment
conditions = NULL : smallint   # Number of conditions 
trials = NULL : smallint       # Number of trials 
matlab = NULL : varchar(100)  # Matlab version used to run the experiment
ptb = NULL: varchar(100)     # Psychtoolbox version used to run the experiment
ns = NULL  : varchar(100)   # Neurosti version used to run the experiment
run = NULL : smallint       # How many times this subject has run this experiment
seq = NULL : smallint       # The sequential recruitment number of this subject for this experiment
%}
%
% BK - April 2022.

classdef Experiment  < dj.Manual
    methods (Access = public)
        function v = file(tbl,root)
            % Return the full filename for the experiments in this table
            data = fetch(tbl*ns.Session,'session_date','file');
            dates = cellfun(@(x) datestr((x),'YYYY/mm/DD'),{data.session_date},'UniformOutput',false); %#ok<DATST>
            v = string(fullfile(dates',{data.file}'));
            if nargin <2
                root =getenv('NS_ROOT');
            end
            v= fullfile(root,v);
        end
        function obj = open(tbl,varargin)
            % Read the experiment files and return an array of Neurostim CIC objects.
            % For standard Neurostim PTB experiments the data are stored as
            % a single c variable in a .mat file. This is loaded using
            % the Matlab buitlin load()
            % By passing a 'fun', (a function that takes a filename and
            % returns an object) the data can be read differently, or
            % preprocessed.
            p=inputParser;
            p.addParameter('fun',@(x)(ns.Experiment.load(x)))
            p.addParameter('mapRoot',{},@iscellstr);
            p.addParameter('perFile',true,@islogical);
            p.parse(varargin{:});
            obj=[];
            if p.Results.perFile
                % Load each file separately and return a vector of objects
                for key=tbl.fetch('file')'
                    if ~isempty(p.Results.mapRoot)
                        thisFile = replace(key.file,p.Results.mapRoot{:});
                    else
                        thisFile = key.file;
                    end
                    obj = [obj; p.Results.fun(thisFile)]; %#ok<AGROW>
                end
            else
                % Load all files at once (only works if the 'fun' passed
                % here can handle that).
                keys = tbl.fetch('file');
                files= {keys.file};
                if ~isempty(p.Results.mapRoot)
                    thisFile = replace(files,p.Results.mapRoot{:});
                else
                    thisFile = files;
                end
                obj = p.Results.fun(thisFile);
            end
        end
        function [out,filename] = get(tbl,varargin)
            % function [out,filename] = get(o,varargin)
            % Retrieve all information on a specific plugin in an experiment
            %
            % INPUT
            % o - A ns.Experiment table with at least 1 row (i.e. one
            %           experiment)
            % plg - A cell array with names of plugins. The {} cell array
            % will retrieve parameters for all plugins. Note that the .cic
            % plugin parameters are always included, whether requested or
            % not.
            % OUTPUT
            % out -  A struct with fields named after the plugins, and each
            % field is another struct that continas all parameter info.
            % All global constants (i.e. those that do
            % not change within an experiment), parameters (a single value per
            % trial that does not change within a trial), and events (which can
            % happen at any time).
            % filename - The file that originally provided these values to the
            % database.
            %
            if ~exists(tbl)
                out = struct([]);
                filename = '';
                return;
            end
            p = inputParser;
            p.addRequired('tbl')
            p.addOptional('plg',{},@(x) (isstring(x) | ischar(x) | iscell(x)))
            p.addParameter('prm',"",@(x) (isstring(x) | ischar(x)))
            p.addParameter('atTrialTime',[],@isnumeric);
            p.parse(tbl,varargin{:});

            if ischar(p.Results.plg)
                plg = {p.Results.plg};
            else
                plg = p.Results.plg;
            end


            ix =1:count(tbl);
            out = cell(numel(ix),1);
            filename = cell(numel(ix),1);
            cntr =0;
            for exptKey=tbl.fetch()'
                cntr = cntr + 1;
                filename{cntr} = fetch1(tbl &exptKey,'file');
                if nargin <2 || isempty(plg)
                    % Get info from all plugins
                    plg  = fetchn( (tbl & exptKey) * ns.Plugin,'plugin_name');
                end
                % Always get cic
                v.cic =  get(ns.PluginParameter  & (ns.Plugin * (tbl & exptKey) & 'plugin_name=''cic''')) ;
                plg = setdiff(plg,{'cic'});
                for pIx = 1:numel(plg)
                    plgName = plg{pIx};
                    % Get the properties for this plugin
                    if isempty(p.Results.prm)
                        % All prms
                        parms =  ns.PluginParameter  & (ns.Plugin * (tbl & exptKey) & ['plugin_name=''' plgName '''']) ;
                    else
                        % One prm
                        parms =  (ns.PluginParameter & ['property_name=''' p.Results.prm '''']) & (ns.Plugin * (tbl & exptKey) & ['plugin_name=''' plgName '''']) ;
                    end
                    if ~exists(parms)
                        continue;
                    end
                    v.(plgName) = get(parms);
                    % Post-process non-CIC plugins if requested
                    if ~isempty(p.Results.atTrialTime)
                        out{cntr} = ns.attrialtime(v.(plgName),p.Results.prm,p.Results.atTrialTime,v.cic);
                    end
                end
                if ~isempty(p.Results.prm) && ~isempty(p.Results.atTrialTime)
                    % Return only the values, not the time/trial; they are
                    % stored in the out cell array already
                else
                    % Return struct with all info
                    out{cntr}=v;
                end
            end

            % Convenience; remove the wrappping cell if it only a single
            % experiment was queried.
            if numel(ix)==1
                out = out{1};
                filename=filename{1};
            end


        end


        function updateWithFileContents(tbl,varargin)
            % function updateWithFileContents(self,oldKey,newOnly)
            % Read neurostim files to fill the database with the
            % information contained in plugins/stimuli. This can be done
            % automatically (by nsAddToDatajoint), or manually to update information
            % INPUT
            % tbl - (A subset of) the ns.Experiment table to update.
            % oldKey - The primary key of the experiment to update (if not
            % specified or empty, all experiments in tbl will be updated).
            % newOnly  - Set to true to update only those experiments that
            % have no information in the database currently. [true]
            p =inputParser;
            p.addRequired('tbl');
            p.addParameter('oldKey',struct([]));
            p.addParameter('newOnly',true,@islogical);
            p.addParameter('root',getenv('NS_ROOT'));
            p.parse(tbl,varargin{:});

            if isempty(p.Results.oldKey)
                % Run all
                for key=tbl.fetch()'
                    try
                        updateWithFileContents(tbl,'oldKey',key,'newOnly',p.Results.newOnly,'root',p.Results.root);
                    catch me
                        fprintf(2,'Failed updating contents from %s\n (%s)',key.file,me.message);
                    end
                end
                return;
            else
                oldKey = p.Results.oldKey;
            end

            % Check if this eperiment already has file data, if newOnly
            % is true, we skip those.
            if p.Results.newOnly && exists(tbl & oldKey) && ~isnan(fetchn(tbl & oldKey,'stimuli'))
                return;
            end

            % Read the file to add details
            oldTuple = fetch(tbl & oldKey,'paradigm','file');

            lastwarn(''); % Reset

            file = fullfile(p.Results.root,strrep(oldTuple.session_date,'-','/'),oldTuple.file);
            % Now add file contents information
            % If the file cannot be read fully (for instance because some
            % classes are not on the current Matlab path) the cic object
            % will be incomplete and some of the code below will fail.
            try
                c  = ns.Experiment.load(file);
            catch
                % Make sure that the experiment key is in the table
                fprintf(2,'Failed to load %s \n',file)
                return;
            end

            if isfield(c,'ptbVersion')
                ptbVersion = c.ptbVersion.version;
            else
                ptbVersion = 'unknown'; % Early versions did not store this
            end
            actualNrTrialsStarted = max([c.prms.trial.log{:}]);
            if actualNrTrialsStarted <1
                % Cannot read some information in a file without trials..
                % Just putting zeros.
                fprintf('Skipping %s - no completed trials\n',file);
                tuple = struct('stimuli',0,'blocks',0,...
                    'conditions',0,'trials',actualNrTrialsStarted,...
                    'matlab',c.matlabVersion,'ptb',ptbVersion,...
                    'ns','#','run',0,'seq',0);
            else
                % Pull the top level information to put in the tbl
                if isempty(c.runNr)
                    runNr = NaN;
                else
                    runNr = c.runNr;
                end
                if isempty(c.seqNr)
                    seqNr = NaN;
                else
                    seqNr = c.seqNr;
                end

                tuple =struct('stimuli',c.nrStimuli,'blocks',c.nrBlocks,...
                    'conditions',c.nrConditions,'trials',actualNrTrialsStarted,...
                    'matlab',c.matlabVersion,'ptb',ptbVersion,'ns','#','run',runNr,'seq',seqNr);
            end

            % Remove current tuple
            if exists(tbl & oldKey)
                del(tbl & oldKey)
            end

            newTuple = mergestruct(oldTuple,tuple);
            insert(tbl,newTuple);

            % Remove the current plugin info and store currently read
            % information.
            if exists(ns.Plugin & oldKey)
                del(ns.Plugin & oldKey)
            end
            if actualNrTrialsStarted>1
                % re-add each plugin (pluginOrder includes stimuli)
                for plg = [c.pluginOrder c]
                    try
                        plgKey = struct('starttime',oldKey.starttime,'session_date',oldKey.session_date,'subject',oldKey.subject);
                        make(ns.Plugin,plgKey,plg);
                    catch
                        fprintf('Failed to add plugin %s information\n',plg.name);
                    end
                end
            end
        end

    end
    methods (Static)
        function o = load(filename)
            % Default method to open a Neurostim data file (a .mat file
            % containing a CIC class object in the variable 'c'
            s  = load(filename,'c');
            o=s.c;
        end
    end
end
