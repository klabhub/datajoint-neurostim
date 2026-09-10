classdef TestDataJointPipeline < matlab.unittest.TestCase
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
    methods (Test)
        function continuousDataIsPopulated(testCase)
            testCase.report('continuousDataIsPopulated: starting');
            key = testCase.experimentKey();
            c = ns.C & key;
            testCase.report('checking ns.C row, channels, time, rate, and signal');
            testCase.verifyEqual(count(c),1); testCase.verifyEqual(c.channels,[1;2]);
            testCase.verifyEqual(c.time,(0:9)'); testCase.verifyEqual(c.samplingRate,1000);
            testCase.verifyEqual(fetch1(ns.CChannel & key & 'channel=1','signal'),single((0:9)'+1000));
            testCase.report('continuousDataIsPopulated: complete');
        end
        function epochsArePopulatedAndAligned(testCase)
            testCase.report('epochsArePopulatedAndAligned: starting');
            key = testCase.experimentKey();
            testCase.report('populating ns.Epoch'); populate(ns.Epoch & key);
            epoch = ns.Epoch & key;
            testCase.report('checking epoch metadata and six EpochChannel rows');
            testCase.verifyEqual(count(epoch),1); testCase.verifyEqual(epoch.time,(-2:2)'/1000);
            testCase.verifyEqual(epoch.samplingRate,1250); testCase.verifyEqual(count(ns.EpochChannel & key),6);
            testCase.report('checking trial 2/channel 1 extracted signal and onset');
            testCase.verifyEqual(fetch1(ns.EpochChannel & key & 'trial=2' & 'channel=1','signal'),double((3:7)'+1000),'AbsTol',1e-10);
            testCase.verifyEqual(fetch1(ns.EpochChannel & key & 'trial=2' & 'channel=1','onset'),5);
            testCase.report('epochsArePopulatedAndAligned: complete');
        end
        function epochKeySourceRequiresDimensionConditions(testCase)
            testCase.report('epochKeySourceRequiresDimensionConditions: starting');
            key = testCase.experimentKey();
            testCase.report('checking Epoch key source and unpopulated Epoch state');
            epochTable = ns.Epoch; testCase.verifyEqual(count(epochTable.keySource & key),1);
            testCase.report('epochKeySourceRequiresDimensionConditions: complete');
        end
        function tepochIsPopulatedWithoutAveraging(testCase)
            testCase.report('tepochIsPopulatedWithoutAveraging: starting');
            key = testCase.experimentKey();
            tepochKey = mergestruct(key,struct('ttag','syntheticTepoch'));
            testCase.report('ensuring ns.Epoch and ns.EpochChannel are populated');
            populate(ns.Epoch & key);

            if count(ns.EpochChannel & key)==0
                testCase.report('recreating epoch part rows after a prior destructive test');
                previousSafeMode = dj.config('safemode');
                restoreSafeMode = onCleanup(@() dj.config('safemode',previousSafeMode));
                dj.config('safemode',false);
                delete(ns.Epoch & key);
                populate(ns.Epoch & key);
            end

            testCase.report('replacing epoch signals with a known one-cycle sinusoid');
            sinusoid = sin(2*pi*(0:4)'/5);
            sourceRows = fetch(ns.EpochChannel & key,'subject','session_date','starttime','ctag','dimension','etag','filename','paradigm','channel','trial');
            for iRow = 1:numel(sourceRows)
                update(ns.EpochChannel & sourceRows(iRow),'signal',sinusoid);
            end
            testCase.report('inserting an FFT TepochParm with averaging disabled');
            insert(ns.TepochParm,struct('ttag','syntheticTepoch','etag','syntheticEpoch', ...
                'fun',struct('fft',{{}}),'window',[-2 2], ...
                'channels',[1 2],'trials',{{1 2 3}},'conditions',{{}},'average',{{}}));

            testCase.report('populating ns.Tepoch without channel/trial averaging');
            populate(ns.Tepoch & tepochKey);
            tepoch = ns.Tepoch & tepochKey;
            testCase.verifyEqual(count(tepoch),2);
            testCase.verifyEqual(sort(string(fetchn(tepoch,'dependent'))),["amplitude";"phase"]);
            testCase.verifyEqual(fetch1(tepoch & 'dependent="amplitude"','independent'),'frequency');
            testCase.verifyEqual(fetch1(tepoch & 'dependent="amplitude"','x'),[0 250 500]);

            testCase.report('checking one transformed row per source trial/channel');
            tc = ns.TepochChannel & tepochKey;
            testCase.verifyEqual(count(tc),12);
            testCase.verifyEqual(sort(tc.channels),[1;2]);
            testCase.verifyEqual(sort(fetchn(tc,'trial')),repelem((1:3)',4));
            testCase.verifyEqual(unique(fetchn(tc,'nrtrials')),1);
            testCase.verifyEqual(unique(fetchn(tc,'nrchannels')),1);
            amplitude = fetch1(tc & 'dependent="amplitude"' & 'trial=2' & 'channel=1','y');
            testCase.verifyEqual(amplitude,[0;sqrt(5);0],'AbsTol',1e-10);
            phase = fetch1(tc & 'dependent="phase"' & 'trial=2' & 'channel=1','y');
            testCase.verifyEqual(phase(2),-pi/2,'AbsTol',1e-10);
            testCase.verifySize(amplitude,[3 1]);
            testCase.report('tepochIsPopulatedWithoutAveraging: complete');
        end
        function chunkedDeleteRemovesPartRowsInBatches(testCase)
            testCase.report('chunkedDeleteRemovesPartRowsInBatches: starting');
            key = testCase.experimentKey();
            testCase.report('ensuring ns.Epoch and ns.EpochChannel are populated');
            populate(ns.Epoch & key);
            if count(ns.Tepoch & key)>0
                testCase.report('removing Tepoch rows before deleting source EpochChannel rows');
                previousSafeMode = dj.config('safemode');
                restoreSafeMode = onCleanup(@() dj.config('safemode',previousSafeMode));
                dj.config('safemode',false);
                delete(ns.Tepoch & key);
            end
            populate(ns.Epoch & key);
            before = count(ns.EpochChannel & key);
            testCase.verifyEqual(before,6);
            previousLimit = getenv('NS_MAXUNCHUNKEDDELETE');
            restoreLimit = onCleanup(@() setenv('NS_MAXUNCHUNKEDDELETE',previousLimit));
            previousSafeMode = dj.config('safemode');
            restoreSafeMode = onCleanup(@() dj.config('safemode',previousSafeMode));
            dj.config('safemode',false);
            setenv('NS_MAXUNCHUNKEDDELETE','1');
            testCase.report('deleting %d EpochChannel rows with batch size 1',before);
            chunkedDelete(ns.EpochChannel & key,1,{},false);
            testCase.verifyEqual(count(ns.EpochChannel & key),0);
            testCase.verifyEqual(count(ns.Epoch & key),1);
            clear restoreLimit
            testCase.report('chunkedDeleteRemovesPartRowsInBatches: complete');
        end
    end
    methods (Access = private)
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
