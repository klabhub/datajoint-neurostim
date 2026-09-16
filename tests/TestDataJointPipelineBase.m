classdef (Abstract) TestDataJointPipelineBase < matlab.unittest.TestCase
    % Integration tests for the ns.C -> ns.Dimension -> ns.Epoch path.
    % These tests are skipped when no DataJoint server is available.
    properties
        projectRoot; dataRoot; databaseName; oldPath; oldFolder; oldRoot
        schemaCreated = false
    end
    properties (Constant, Access = private)
        subject = 'synthetic1'; sessionDate = '2024-02-29'; startTime = '09:00:00'
    end
    methods (TestClassSetup)
        function createDatabase(testCase)
            testCase.report('class setup: starting');
            testCase.oldPath = path; testCase.oldFolder = pwd;
            testCase.oldRoot = getenv('NS_ROOT');
            testCase.dataRoot = string(tempname); testCase.projectRoot = string(tempname);
            mkdir(testCase.dataRoot); mkdir(testCase.projectRoot);
            testCase.report('temporary data root: %s',testCase.dataRoot);
            repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
            addpath(repositoryRoot,fullfile(repositoryRoot,'datajoint-matlab'), ...
                fullfile(repositoryRoot,'tools','mym',mexext));
            try
                testCase.report('checking DataJoint connection');
                dj.conn;
            catch ME
                testCase.cleanup();
                testCase.assumeTrue(false,['DataJoint integration tests require a connection: ' ME.message]);
                return
            end
            testCase.databaseName = "dj_test_" + lower(string(char(java.util.UUID.randomUUID.toString))).replace('-','');
            testCase.schemaCreated = true; % Arm cleanup before any database operation.
            testCase.report('creating isolated database/schema: %s',testCase.databaseName);
            try
                nsInitializeDataJoint(testCase.projectRoot,testCase.databaseName,"ns",dataRoot=testCase.dataRoot);
            catch ME
                testCase.report('schema setup failed; cleaning up partial setup');
                testCase.cleanup();
                testCase.assumeTrue(false,['Cannot create isolated DataJoint test database: ' ME.message]);
                return
            end
            addpath(testCase.projectRoot); clear ns.getSchema; setenv('NS_ROOT',testCase.dataRoot);
            testCase.report('seeding synthetic subject/session/experiment metadata');
            testCase.seedFixture();
            testCase.report('populating prerequisites for EpochParm validation');
            key = testCase.experimentKey();
            populate(ns.File & key); populate(ns.C & key); populate(ns.Dimension & key);
            testCase.insertEpochParm();
            testCase.report('class setup: complete');
        end
    end
    methods (TestClassTeardown)
        function removeDatabase(testCase)
            testCase.report('class teardown: starting');
            testCase.cleanup();
            testCase.report('class teardown: complete');
        end
    end
    methods (Access = protected)
        function seedFixture(testCase)
            dayFolder = fullfile(testCase.dataRoot,'2024','02','29'); mkdir(dayFolder);
            fileName = 'synthetic.synthetic.090000.mat';
            fid = fopen(fullfile(dayFolder,fileName),'w'); fclose(fid);
            key = testCase.experimentKey();
            testCase.report('inserting Subject, Session, and Experiment');
            insert(ns.Subject,struct('subject',testCase.subject));
            insert(ns.Session,struct('subject',testCase.subject,'session_date',testCase.sessionDate));
            insert(ns.Experiment,mergestruct(key,struct('paradigm','synthetic','file',fileName,'nrtrials',3)));
            testCase.report('inserting synthetic CParm and plugin parameters');
            insert(ns.CParm,struct('ctag','synthetic','fun','testsupport.syntheticCReader','extension','.mat', ...
                'description','test reader','parms',struct('nSamples',10,'startTime',0,'stopTime',9,'channels',[1 2]),'include',fileName));
            insert(ns.Plugin,mergestruct(key,struct('plugin_name','cic')));
            insert(ns.Plugin,mergestruct(key,struct('plugin_name','synthetic')));
            cicKey = mergestruct(key,struct('plugin_name','cic'));
            addNew(ns.PluginParameter,cicKey,'firstframe',[0 1 9],'Event',[0 0 0],[1 2 3],[0 1 9]);
            pluginKey = mergestruct(key,struct('plugin_name','synthetic'));
            addNew(ns.PluginParameter,pluginKey,'startTime',[2 5 8],'Event',[2 4 1],[1 2 3],[2 5 8]);
            addNew(ns.PluginParameter,pluginKey,'condition',[1 2 1],'Parameter',[0 0 0],[1 2 3],[2 5 8]);
            testCase.report('inserting DimensionParm and EpochParm');
            insert(ns.DimensionParm,struct('dimension','condition','paradigm','synthetic', ...
                'parms',struct('plg','synthetic','prm','condition','atTrialTime',0)));
        end
        function cleanup(testCase)
            if testCase.schemaCreated && ~isempty(testCase.databaseName)
                try
                    testCase.report('dropping database if it exists: %s',testCase.databaseName);
                    dj.config('safemode',false);
                    query(dj.conn,sprintf('DROP DATABASE IF EXISTS `%s`',testCase.databaseName));
                catch ME
                    warning('TestDataJointPipeline:TeardownFailed','Could not remove test database: %s',ME.message);
                end
            end
            setenv('NS_ROOT',testCase.oldRoot);
            if ~isempty(testCase.oldPath), path(testCase.oldPath); end
            if ~isempty(testCase.oldFolder) && isfolder(testCase.oldFolder), cd(testCase.oldFolder); end
            if ~isempty(testCase.dataRoot) && isfolder(testCase.dataRoot), testCase.report('removing temporary data root'); rmdir(testCase.dataRoot,'s'); end
            if ~isempty(testCase.projectRoot) && isfolder(testCase.projectRoot), testCase.report('removing temporary project root'); rmdir(testCase.projectRoot,'s'); end
            testCase.schemaCreated = false;
        end
        function report(~,message,varargin)
            fprintf('[TestDataJointPipeline] %s\n',sprintf(message,varargin{:}));
        end
        function insertEpochParm(~)
            insert(ns.EpochParm,struct('etag','syntheticEpoch','ctag','synthetic','dimension','condition','window',[-2 2], ...
                'align',struct('plugin','synthetic','event','startTime'),'prepparms',struct('enable',false), ...
                'artparms',struct('enable',false),'plgparms',struct('enable',false)));
        end
        function key = experimentKey(testCase)
            key = struct('subject',testCase.subject,'session_date',testCase.sessionDate,'starttime',testCase.startTime);
        end
    end
end
