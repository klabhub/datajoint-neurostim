classdef TestDataJointPipelineFixture < matlab.unittest.fixtures.Fixture
    % Shared DataJoint database and synthetic pipeline fixture.
    properties (SetAccess = private)
        projectRoot
        dataRoot
        databaseName
    end
    properties (Access = private)
        oldPath
        oldFolder
        oldRoot
        schemaCreated = false
    end
    methods
        function setup(fixture)
            fixture.oldPath = path;
            fixture.oldFolder = pwd;
            fixture.oldRoot = getenv('NS_ROOT');
            fixture.dataRoot = string(tempname);
            fixture.projectRoot = string(tempname);
            mkdir(fixture.dataRoot);
            mkdir(fixture.projectRoot);
            fixture.report('class setup: starting');
            fixture.report('temporary data root: %s',fixture.dataRoot);
            repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
            addpath(repositoryRoot,fullfile(repositoryRoot,'datajoint-matlab'), ...
                fullfile(repositoryRoot,'tools','mym',mexext));
            addpath(fileparts(mfilename('fullpath')));
            try
                fixture.report('checking DataJoint connection');
                dj.conn;
            catch ME
                fixture.cleanup();
                rethrow(ME);
            end
            fixture.databaseName = "dj_test_" + ...
                lower(string(char(java.util.UUID.randomUUID.toString))).replace('-','');
            fixture.schemaCreated = true;
            fixture.report('creating isolated database/schema: %s',fixture.databaseName);
            try
                nsInitializeDataJoint(fixture.projectRoot,fixture.databaseName,"ns", ...
                    dataRoot=fixture.dataRoot);
            catch ME
                fixture.report('schema setup failed; cleaning up partial setup');
                fixture.cleanup();
                rethrow(ME);
            end
            addpath(fixture.projectRoot);
            clear ns.getSchema;
            setenv('NS_ROOT',fixture.dataRoot);
            fixture.seedFixture();
            key = fixture.experimentKey();
            fixture.report('populating prerequisites for EpochParm validation');
            populate(ns.File & key);
            populate(ns.C & key);
            populate(ns.Dimension & key);
            fixture.insertEpochParm();
            fixture.report('class setup: complete');
        end
        function teardown(fixture)
            fixture.report('class teardown: starting');
            fixture.cleanup();
            fixture.report('class teardown: complete');
        end
    end
    methods (Access = private)
        function seedFixture(fixture)
            dayFolder = fullfile(fixture.dataRoot,'2024','02','29');
            mkdir(dayFolder);
            fileName = 'joe.synthetic.090000.mat';
            fid = fopen(fullfile(dayFolder,fileName),'w');
            fclose(fid);
            key = fixture.experimentKey();
            fixture.report('inserting Subject, Session, and Experiment');
            insert(ns.Subject,struct('subject','joe'));
            insert(ns.Session,struct('subject','joe','session_date','2024-02-29'));
            insert(ns.Experiment,mergestruct(key,struct( ...
                'paradigm','synthetic','file',fileName,'nrtrials',3)));
            fixture.report('inserting synthetic CParm and plugin parameters');
            insert(ns.CParm,struct('ctag','synthetic', ...
                'fun','testsupport.syntheticCReader','extension','.mat', ...
                'description','test reader','parms',struct( ...
                'nSamples',10,'startTime',0,'stopTime',9,'channels',[1 2]), ...
                'include',fileName));
            insert(ns.Plugin,mergestruct(key,struct('plugin_name','cic')));
            insert(ns.Plugin,mergestruct(key,struct('plugin_name','synthetic')));
            cicKey = mergestruct(key,struct('plugin_name','cic'));
            addNew(ns.PluginParameter,cicKey,'firstframe',[0 1 9], ...
                'Event',[0 0 0],[1 2 3],[0 1 9]);
            pluginKey = mergestruct(key,struct('plugin_name','synthetic'));
            addNew(ns.PluginParameter,pluginKey,'startTime',[2 5 8], ...
                'Event',[2 4 1],[1 2 3],[2 5 8]);
            addNew(ns.PluginParameter,pluginKey,'condition',[1 2 1], ...
                'Parameter',[0 0 0],[1 2 3],[2 5 8]);
            fixture.report('inserting DimensionParm and EpochParm');
            insert(ns.DimensionParm,struct('dimension','condition', ...
                'paradigm','synthetic','parms',struct( ...
                'plg','synthetic','prm','condition','atTrialTime',0)));
        end
        function insertEpochParm(~)
            insert(ns.EpochParm,struct('etag','syntheticEpoch', ...
                'ctag','synthetic','dimension','condition','window',[-2 2], ...
                'align',struct('plugin','synthetic','event','startTime'), ...
                'prepparms',struct('enable',false), ...
                'artparms',struct('enable',false), ...
                'plgparms',struct('enable',false)));
        end
        function key = experimentKey(~)
            key = struct('subject','joe', ...
                'session_date','2024-02-29', ...
                'starttime','09:00:00');
        end
        function cleanup(fixture)
            if fixture.schemaCreated && ~isempty(fixture.databaseName)
                try
                    fixture.report('dropping database if it exists: %s',fixture.databaseName);
                    dj.config('safemode',false);
                    query(dj.conn,sprintf('DROP DATABASE IF EXISTS %s',fixture.databaseName));
                catch ME
                    warning('TestDataJointPipeline:TeardownFailed', ...
                        'Could not remove test database: %s',ME.message);
                end
            end
            setenv('NS_ROOT',fixture.oldRoot);
            if ~isempty(fixture.oldPath), path(fixture.oldPath); end
            if ~isempty(fixture.oldFolder) && isfolder(fixture.oldFolder)
                cd(fixture.oldFolder);
            end
            if ~isempty(fixture.dataRoot) && isfolder(fixture.dataRoot)
                fixture.report('removing temporary data root');
                rmdir(fixture.dataRoot,'s');
            end
            if ~isempty(fixture.projectRoot) && isfolder(fixture.projectRoot)
                fixture.report('removing temporary project root');
                rmdir(fixture.projectRoot,'s');
            end
            fixture.schemaCreated = false;
        end
        function report(~,message,varargin)
            fprintf('[TestDataJointPipeline] %s\n',sprintf(message,varargin{:}));
        end
    end
end
