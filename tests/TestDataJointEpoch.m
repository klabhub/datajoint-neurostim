classdef TestDataJointEpoch < TestDataJointPipelineBase
    methods (Test)
        function windowMaskingRetainsTrialsAndStoresIntervals(testCase)
            key = testCase.experimentKey();
            etag = 'syntheticMaskedEpoch';
            options = struct('enable',true,'action',"nan",'window',2, ...
                'step',1,'padding',0,'amplitude_threshold_peak',0.5, ...
                'variance_z_threshold',1e9,'hf_z_threshold',1e9, ...
                'correlation_z_threshold',1e9);
            insert(ns.EpochParm,struct('etag',etag,'ctag','synthetic', ...
                'dimension','condition','window',[-2 2], ...
                'align',struct('plugin','synthetic','event','startTime'), ...
                'artparms',options));
            testCase.addTeardown(@() TestDataJointEpoch.removeTestEpoch(etag));
            key.etag = etag;

            populate(ns.Epoch & key);

            testCase.verifyEqual(count(ns.Epoch & key),1);
            testCase.verifyEqual(count(ns.EpochChannel & key),6);
            signals = fetchn(ns.EpochChannel & key,'signal');
            testCase.verifyTrue(all(cellfun(@(x) all(isnan(x)),signals)));
            art = fetch1(ns.Epoch & key,'art');
            testCase.verifyEmpty(art.badBy.all);
            testCase.verifyEqual(art.segments.maskedFraction,ones(3,1));
            testCase.verifyEqual(sort(art.segments.intervals.trial),(1:3)');
            testCase.verifyTrue(isfield(art.segments.intervals,'start'));
            summary = attrition(ns.Epoch & key);
            testCase.verifyEqual(summary.badByArtifact,0);
            testCase.verifyEqual(summary.artifactMaskedFraction,1);
        end

        function epochsArePopulatedAndAligned(testCase)
            testCase.report('epochsArePopulatedAndAligned: starting');
            key = testCase.experimentKey();
            testCase.report('populating ns.Epoch'); populate(ns.Epoch & key);
            epoch = ns.Epoch & key;
            testCase.report('checking epoch metadata and six EpochChannel rows');
            testCase.verifyEqual(count(epoch),1);
            epochTime = fetch1(epoch,'time');
            testCase.verifyEqual(linspace(epochTime(1),epochTime(2),epochTime(3))',(-2:2)'/1000);
            c = ns.C & key;
            testCase.verifyEqual(c.samplingRate,1000);
            testCase.verifyEqual(count(ns.EpochChannel & key),6);
            signals = fetchn(ns.EpochChannel & key,'signal');
            testCase.verifyTrue(all(cellfun(@isrow,signals)));
            testCase.report('checking trial 2/channel 1 extracted signal and onset');
            testCase.verifyEqual(fetch1(ns.EpochChannel & key & 'trial=2' & 'channel=1','signal')',double((3:7)'+1000),'AbsTol',1e-10);
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
    end
    methods (Static, Access=private)
        function removeTestEpoch(etag)
            oldSafeMode = dj.config('safemode');
            restore = onCleanup(@() dj.config('safemode',oldSafeMode));
            dj.config('safemode',false);
            del(ns.EpochParm & struct('etag',etag));
        end
    end
end
