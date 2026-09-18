classdef TestDataJointEpoch < TestDataJointPipelineBase
    methods (Test)
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
end
