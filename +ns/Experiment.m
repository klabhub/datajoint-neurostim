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
nrtrials = NULL : smallint       # Number of trials 
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
        function [ana,notAna] = analyze(tbl,pv)
            % Return the subtable that is marked to be analyzed in the meta
            % data. The optional second output argument is the table of
            % experiments that will not be analyzed.
            % Use strict = true to select only those experiments in which
            % there is a "1"  entry for the meta value 'analyze', use
            % strict = false to include experiments as long as they do not have
            % an analyze meta field that is "0"  (including undefined analyze meta fields).
            arguments
                tbl (1,1) ns.Experiment
                pv.strict (1,1) logical = false; % Set to true to require an explicit analyze =1 setting
            end
            exptWithMetaAnalyze = tbl*ns.ExperimentMeta & proj(ns.ExperimentMeta & 'meta_name="analyze"');

            if pv.strict
                % Only use the ones explicitly set to "1" (default is
                % analyze = false)
                ana  =  tbl & proj(exptWithMetaAnalyze & 'meta_value ="1"') ;
            else
                % Use the ones NOT set to 0 and the ones where analyze is
                % not defined (i.e. default assumption is analyze=true)
                exptWithout = tbl- proj(exptWithMetaAnalyze);
                ana  =  tbl & (proj(exptWithMetaAnalyze & ('NOT meta_value ="0"')) | proj(exptWithout));
            end
            if nargout>1
                notAna = tbl  -proj(ana);
            end

        end
        function v = folder(tbl,root)
            % Return the full path to the experiments in this table
            if nargin <2
                root =getenv('NS_ROOT');
            end
            v = fileparts(file(tbl,root));
        end
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
                    ff = file(ns.Experiment & key);
                    if ~isempty(p.Results.mapRoot)
                        thisFile = replace(ff,p.Results.mapRoot{:});
                    else
                        thisFile = ff;
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

        function nwbRoot = nwbExport(expt,pv)
            % Create Neurodata without Borders objects for each experiment in a
            % ns.Experiment table.
            %
            % folder = Folder where the data files will be saved. By default this is
            % empty, which will create and return the NwbFile objects, but not save
            % them.
            % tz  - TimeZone to use for session_Start_time and
            % timestamps_reference_time ["local"]
            % general -  A struct with fields named after the general_ properties of the
            % NwbFile object that use free form text to describe meta data.
            % The allowable names depend on the NWM file schema.  By default this is empty.
            %
            %
            % EXAMPLE
            %  general = struct('lab',"KLab",'Insitution',"Rutgers University - Newark")
            %
            % The columns of the ns.Subject table are always included, to include ns.SubjectMeta
            % data as well, provide a dictionary that maps the names of the meta data
            % to names in the NWB subject schema.
            %
            % subjectMeta =dictionary('strain','strain;')
            %
            % Datajoint Tables (classes) with a nwb member function will be called
            % automatically to add their contribution to the nwb object.  The nwb
            % member takes the table as its first input, the root of the nwb object
            % as its second, adnd the pv struct from this function as its third.
            %  Prototype: nwb(dj.Relvar,NwbFile,struct)
            % see ns.Subject for an example
            %
            % createNwb(ns.Experiment ,general=general, subjectMeta= subjectMeta,
            %                               folder = "c:/temp");
            % BK - 2023,2024
            arguments
                expt (1,1) ns.Experiment
                pv.packages (1,:)= ""
                pv.folder (1,1) = ""
                pv.tz (1,1) string = "local"
                pv.general (1,1) struct = struct();
                pv.subjectMeta (1,1) dictionary  = dictionary(string([]),string([]));
            end

            assert(~isempty(which('NwbFile')),'This function depends on the matnwb package. Install it from github and add it to the Matlab path');

            classesWithNwb =nwbFind("ns.Subject",pv.packages);

            for e = fetch(expt,'*')'
                % Loop over experiments (NWB refers to this as a "session", NS uses session to refer to all experiments for
                % a subject on a given day)
                uniqueExperimentName = sprintf('%s_%s_%s',e.subject,e.session_date,e.starttime);
                % The start time of the experiment:  (note that this is not
                % exactly the same as the timestapmps_reference_time;the
                % former is based on a call to now, while the latter is a
                % call to GetSecs that is executed a few lines later in the
                % constructor of cic.
                startTime = datetime([e.session_date 'T' e.starttime],'TimeZone',pv.tz);

                nwbRoot = NwbFile(...
                    'session_description',sprintf('%s',e.paradigm ),...
                    'identifier',char(java.util.UUID.randomUUID().toString()),...
                    'session_start_time', startTime,...
                    'timestamps_reference_time',startTime, ...
                    'general_session_id', uniqueExperimentName);



                %% General properties (specified in the call to this function)
                fn = fieldnames(pv.general);
                for i=1:numel(fn)
                    prop = "general_" + fn{i};
                    try
                        nwbRoot.(prop) = pv.general.(fn{i});
                    catch
                        fprintf(2,"Property %s does not exist in the NWB schema. Ignored. \n", prop );
                    end
                end


                %% Export all tables that have nwb functionality
                for cls=classesWithNwb
                    tbl = feval(cls) & e;
                    nwb(tbl,nwbRoot,pv);
                end


                %% Export to file
                if pv.folder~=""
                    assert(exist(pv.folder,"dir"),"Folder %s does not exist.",pv.folder);
                    fname = fullfile(pv.folder,[strrep(uniqueExperimentName,':','') '.nwb']);
                    nwbExport(nwbRoot,fname);
                end

            end
        end


        function [out,filename] = get(tbl,plg,pv)
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
            % 'prm' -  Specify a single parameter to retrieve
            % 'what' - What is it about the single parameter you want to
            % retrieve.
            %       'data' : the value
            %       'trialtime' : the trial time at which the parameter was
            %       set (0= firstFrame in the trial)
            %       'clocktime' : the time on the neurostim clock when the
            %       parameter was set
            %       'trial' : the trial in which the parameter was set.
            % 'atTrialTime' : Retrieve the value (data) of a parameter at a
            % certain time in each trial. Use this for parameters that
            % define a conditio (e.g., orientation) to get the values at
            % the start of the trial.  (atTrialTime overrules the 'what'
            % parameter and always returns the data).
            %
            %
            % OUTPUT
            % out -  A struct with fields named after the plugins, and each
            % field is another struct that contains all parameter info.
            % This includes all global constants (i.e. those that do
            % not change within an experiment), parameters (a single value per
            % trial that does not change within a trial), and events (which can
            % happen at any time).
            %
            % When 'prm' is used, the output is a vector of values.
            %
            % filename - The file that originally provided these values to the
            % database.
            %
            arguments
                tbl (1,1) ns.Experiment {mustHaveRows}
                plg {mustBeText} =  {''}
                pv.prm {mustBeText} = ''
                pv.what {mustBeNonzeroLengthText,mustBeMember(pv.what,["data" "trialtime" "trial" "clocktime"])} = "data"
                pv.atTrialTime (1,1) double = NaN
                pv.trial (1,:) double = []
            end
            if ischar(plg)
                plg = {plg};
            end

            nrPlugins = sum(cellfun(@(x) (strlength(x)>0),plg,'uni',true));

            ix =1:count(tbl);
            out = cell(numel(ix),1);
            filename = cell(numel(ix),1);
            exptCntr =0;
            for exptKey=tbl.fetch('paradigm')'
                exptCntr = exptCntr + 1;
                filename{exptCntr} = fetch1(tbl &exptKey,'file');
                if nrPlugins==0
                    % Get info from all plugins
                    plg  = fetchn( (tbl & exptKey) * ns.Plugin,'plugin_name');
                    nrPlugins = numel(plg);
                end
                % Always get ***all*** prms from cic
                v.cic =  get(ns.PluginParameter  & (ns.Plugin * (tbl & exptKey) & 'plugin_name=''cic''')) ;
                for pIx = 1:nrPlugins
                    plgName = plg{pIx};
                    if strcmpi(plgName,'cic');continue;end % skip
                    % Get the properties for this plugin
                    if strlength(pv.prm)==0
                        % All prms
                        parms =  ns.PluginParameter  & (ns.Plugin * (tbl & exptKey) & ['plugin_name=''' plgName '''']) ;
                    else
                        % One prm
                        parms =  (ns.PluginParameter & "property_name='"+string(pv.prm) + "'") & (ns.Plugin * (tbl & exptKey) & ['plugin_name=''' plgName '''']) ;
                    end
                    if ~exists(parms)
                        continue;
                    end
                    v.(plgName) = get(parms);
                end

                notFound = ~isfield(v,plg);
                if any(notFound)
                    fprintf('This experiment (%s:%s/%s@%s) did not use the %s plugin(s) or the plugin does not have the reuqested parameter.\n',exptKey.paradigm,exptKey.subject,exptKey.session_date,exptKey.starttime, strjoin(plg(notFound),'/'));
                    out{exptCntr} = [];
                    continue; % Next experiment
                end


                if strlength(pv.prm) ~=0  %prm was specified
                    % Single prm from a specified plugin,  at a specific time
                    out{exptCntr} = ns.attrialtime(v.(plg{1}),pv.prm,pv.atTrialTime,v.cic,pv.what,pv.trial);
                else
                    % Everything
                    % Return struct with all info
                    out{exptCntr}=v;
                end
            end

            % Convenience; remove the wrappping cell if it only a single
            % experiment was queried.
            if numel(ix)==1
                out = out{1}; % Single column
                filename=filename{1};
            end

            if strlength(pv.prm)~=0
                % Single parm, join values as columns
                if iscellstr(out)
                    out = out(:)'; %trials as columns
                elseif iscell(out)
                    if all(cellfun(@numel,out)==numel(out{1}))
                        out =cat(2,out{:}); %trials as columns
                    else
                        out = out'; % For trials as columns
                        return
                    end
                end
                if iscell(out) && ~iscellstr(out)
                    out = cell2mat(out);
                end
            end


        end


        function updateWithFileContents(tbl,cic,pv)
            % function updateWithFileContents(self)
            % Read neurostim files to fill the database with the
            % information contained in plugins/stimuli. This can be done
            % automatically (by nsScan/nsAddToDataJoint), or manually
            % to update information for experiments that were initially
            % added without reading file contents.
            %
            % INPUT
            % tbl - (A subset of) the ns.Experiment table to update.
            % cic - A vector of cic objects already loaded (by nsScan) and
            %           corresponding to the rows of the tbl.
            % newOnly  - Set to true to update only those experiments that
            % have no information in the database currently. [true]
            arguments
                tbl ns.Experiment
                cic (1,:) = []
                pv.newOnly (1,1) logical = true
                pv.pedantic (1,1) logical = false
            end
            tic;
            fprintf('Updating ns.Experiment with file contents from %d experiments...\n',count(tbl))
            % Run all
            keyCntr = 0;
            pkey = tbl.primaryKey;

            for key=tbl.fetch('file')'
                keyCntr=keyCntr+1;
                if isempty(cic)
                    % Read from files
                    % Check if this eperiment already has file data, if newOnly
                    % is true, we skip those.
                    if pv.newOnly && ~isnan(fetchn(tbl & key,'stimuli'))
                        continue;
                    end
                    [thisTpl,thisC] = ns.Experiment.readCicContents(key);
                else
                    % Cic was passed check that it matches, then add
                    % contents. The order of the tbl is not guaranteed, so make sure
                    % to matchup with the correct cic.
                    cicUID = string(datetime({cic.date},'InputFormat','dd MMM yyyy','Format','yyyy-MM-dd'))+string({cic.file})+".mat"; 
                    stay = cicUID==string([key.session_date key.file]);
                    thisC = cic(stay);
                    thisTpl = ns.Experiment.tplFromCic(thisC);
                    thisTpl =mergestruct(key,thisTpl); % Errors if thisC does not belong to this key.
                end

                % Remove current Experiment tuple
                if pv.pedantic
                    if exists(tbl&key)
                        del(tbl & key);
                    end
                    insert(tbl,thisTpl);
                else
                    %Update each field. Potential for referential integrity
                    %loss (not in practice, though).
                    fieldsToUpdate = setdiff(fieldnames(thisTpl),pkey)';
                    % The update could change the tbl (if defined as a
                    % query that uses the fields to be updated), so use
                    % the full ns.Experiment  with a key restriction
                    % instead
                    for fld = fieldsToUpdate
                        update(ns.Experiment &key,fld{1},thisTpl.(fld{1}))
                    end
                end


                % Remove the current plugin info and store currently read
                % information.
                if exists(ns.Plugin & key)
                    del(ns.Plugin & key);
                end
                if max([thisC.prms.trial.log{:}])>0
                    % re-add each plugin (pluginOrder includes stimuli)
                    plgsToAdd= [thisC.pluginOrder thisC];
                    plgKey = struct('starttime',thisTpl.starttime,'session_date',thisTpl.session_date,'subject',thisTpl.subject);
                    for plg = plgsToAdd
                        try
                            make(ns.Plugin,plgKey,plg);
                        catch me
                            fprintf(2,'Failed to add plugin %s information\n (%s)',plg.name,me.message);
                        end
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

        function [tpl,c] = readCicContents(tpl,pv)
            % This is called by nsScan to read the contents of a neurostim
            % output file and extract some meta data to store in DataJoint
            % or show in nsMeta. This is a static so that it can be used
            % with nsMeta without datajoint access.
            arguments
                tpl
                pv.root {mustBeText} = getenv('NS_ROOT');
            end

            file = fullfile(pv.root,strrep(tpl.session_date,'-','/'),tpl.file);
            % If the file cannot be read fully (for instance because some
            % classes are not on the current Matlab path) the cic object
            % will be incomplete and some of the code below will fail.

            % Turn off this warning as it is more or less expected that
            % some will fail, but we're ignoring them anyway
            warnstate= warning('query');
            warning('off','MATLAB:load:classNotFound');
            warning('off','MATLAB:class:DefaultObjectSubstitution');
            try
                c  = ns.Experiment.load(file);
            catch
                % Make sure that the experiment key is in the table
                fprintf(2,'Failed to load %s . Is NS_ROOT (%s) correct? \n',file,getenv('NS_ROOT'))
                return;
            end
            warning(warnstate)

            tpl = ns.Experiment.tplFromCic(c);
        end

        function tpl = tplFromCic(c)
            % Construct a tpl from a CIC object
            arguments
                c (1,1) neurostim.cic
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
                fprintf('Skipping %s - no completed trials\n',c.file);
                tpl = struct('stimuli',0,'blocks',0,...
                    'conditions',0,'nrtrials',actualNrTrialsStarted,...
                    'matlab',c.matlabVersion,'ptb',ptbVersion,...
                    'ns','#','run',0,'seq',0,'paradigm',c.paradigm,'file',[c.file '.mat']);
            else
                % Pull the top level information to put in the tbl
                if isempty(c.runNr)
                    runNr =0;  % early cic did not define this. Default to 0 ( nan or [] causes mysql trouble)
                else
                    runNr = c.runNr;
                end
                if isempty(c.seqNr)
                    seqNr = 0;% early cic did not define this. Default to 0 ( nan or [] causes mysql trouble)
                else
                    seqNr = c.seqNr;
                end

                tpl =struct('stimuli',c.nrStimuli,'blocks',c.nrBlocks,...
                    'conditions',c.nrConditions,'nrtrials',actualNrTrialsStarted,...
                    'matlab',c.matlabVersion,'ptb',ptbVersion,'ns','#','run',runNr,'seq',seqNr,'paradigm',c.paradigm,'file',[c.file '.mat']);
            end
            key =struct('starttime',c.startTimeStr,'session_date',char(datetime(c.date,'InputFormat','dd MMM yyyy','Format','yyyy-MM-dd')),'subject',c.subject);
            tpl  =mergestruct(key,tpl);
        end

    end
end
