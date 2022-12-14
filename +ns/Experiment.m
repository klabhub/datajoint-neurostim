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
        function v = file(tbl)
            v = string(fetchn(tbl ,'file'));
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


        function updateWithFileContents(tbl,oldKey,newOnly)
            % function updateWithFileContents(self,oldKey,newOnly)
            % Read neurostim files to fill the database with the
            % information contained in plugins/stimuli. This can be done
            % automatically (by ns.scan), or manually to update information
            % INPUT
            % tbl - (A subset of) the ns.Experiment table to update.
            % oldKey - The primary key of the experiment to update (if not
            % specified or empty, all experiments in tbl will be updated).
            % newOnly  - Set to true to update only those experiments that
            % have no information in the database currently. [true]

            if nargin <3
                newOnly = true;
            end
            if nargin<2 || isempty(oldKey)
                % Run all
                for key=tbl.fetch()'
                    updateWithFileContents(tbl,key,newOnly);
                end
                return;
            end
            % Check if this eperiment already has file data, if newOnly
            % is true, we skip those.
            if newOnly && exists(tbl & oldKey) && ~isnan(fetchn(tbl & oldKey,'stimuli'))
                return;
            end

            % Read the file to add details
            if isfield(oldKey,'file')
                file = oldKey.file;
            else
                file = fetch1(tbl & oldKey, 'file');
            end
            lastwarn(''); % Reset
            c  = ns.Experiment.load(file);
            if contains(lastwarn,'Cannot load an object of class')
                fprintf(2,'Please add the path to all Neurostim classes to your path.\n')
                fprintf(2,'Skipping %s\n',file)
            else
                actualNrTrialsStarted = max([c.prms.trial.log{:}]);
                if actualNrTrialsStarted <1
                    % Cannot read some information in a file without trials..
                    % Just putting zeros.
                    fprintf('Skipping %s - no completed trials\n',file);
                    tuple = {'stimuli',0,'blocks',0,...
                        'conditions',0,'trials',actualNrTrialsStarted,...
                        'matlab',c.matlabVersion,'ptb',c.ptbVersion.version,...
                        'ns','somehash','run',0,'seq',0};
                else
                    % Pull the top level information to put in the tbl
                    tuple = {'stimuli',c.nrStimuli,'blocks',c.nrBlocks,...
                        'conditions',c.nrConditions,'trials',actualNrTrialsStarted,...
                        'matlab',c.matlabVersion,'ptb',c.ptbVersion.version,...
                        'ns','somehash','run',c.runNr,'seq',c.seqNr};
                end
                if ~exists(tbl & oldKey)
                    insert(tbl,oldKey)
                end
                % Then update
                for i=1:2:numel(tuple)
                    update(tbl & oldKey,tuple{i:i+1}); %Updated version
                end
                % Remove the current plugin info.
                if exists(ns.Plugin & oldKey)
                    del(ns.Plugin & oldKey)
                end
                if actualNrTrialsStarted>1
                    % re-add each plugin (pluginOrder includes stimuli)
                    for plg = [c.pluginOrder c]
                        plgKey = struct('starttime',oldKey.starttime,'session_date',oldKey.session_date,'subject',oldKey.subject);
                        make(ns.Plugin,plgKey,plg);
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
