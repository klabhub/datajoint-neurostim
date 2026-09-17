classdef TestDataJointTepoch < TestDataJointPipelineBase
    methods (Test)
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
                'fun',struct('fft',struct('n',5)), ...
                'parms',struct('timeWindow',[-2 2],'channel',[1 2],'trial',[1 2 3], ...
                'average',string.empty)));

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
            amplitude = fetch1(tc & 'dependent="amplitude"' & 'trial=2' & 'channel=1','signal');
            testCase.verifyEqual(amplitude,[0;sqrt(5);0],'AbsTol',1e-10);
            phase = fetch1(tc & 'dependent="phase"' & 'trial=2' & 'channel=1','signal');
            testCase.verifyEqual(phase(2),-pi/2,'AbsTol',1e-10);
            testCase.verifySize(amplitude,[3 1]);
            testCase.report('tepochIsPopulatedWithoutAveraging: complete');
        end
        function tepochAveragesAcrossTrialsAndChannels(testCase)
            key = testCase.experimentKey();
            tepochKey = mergestruct(key,struct('ttag','averagedTepoch'));
            populate(ns.Epoch & key);
            testCase.setKnownEpochSignals(key);
            insert(ns.TepochParm,struct('ttag','averagedTepoch','etag','syntheticEpoch', ...
                'fun',struct('msten',struct([])), ...
                'parms',struct('timeWindow',[-2 2],'channel',[1 2],'trial',[1 2 3], ...
                'average',["trial" "channel"])));
            populate(ns.Tepoch & tepochKey);
            tc = ns.TepochChannel & tepochKey;
            actualMeans = fetchn(tc & struct('dependent','mean'),'signal');
            normalizedActualMeans = cellfun(@(signal) reshape(signal,1,[]), ...
                actualMeans,UniformOutput=false);
            actualMeans = sortrows(cat(1,normalizedActualMeans{:}),1);
            testCase.verifyEqual(count(tc),6);
            testCase.verifyEqual(unique(fetchn(tc,'channel')),0);
            testCase.verifyEqual(unique(fetchn(tc,'trial')),0);
            testCase.verifyEqual(unique(fetchn(tc,'nrtrials')),[1;2]);
            testCase.verifyEqual(unique(fetchn(tc,'nrchannels')),2);
            expectedMeans = [190+(0:4);200+(0:4)];
            testCase.verifyEqual(actualMeans,expectedMeans,'AbsTol',1e-10);
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
    methods (Access=private)
        function setKnownEpochSignals(~,key)
            sourceRows = fetch(ns.EpochChannel & key,'subject','session_date','starttime', ...
                'ctag','dimension','etag','filename','paradigm','channel','trial');
            for iRow = 1:numel(sourceRows)
                signal = 100*sourceRows(iRow).channel + 10*sourceRows(iRow).trial^2 + (0:4)';
                update(ns.EpochChannel & sourceRows(iRow),'signal',signal);
            end
        end
    end
end
