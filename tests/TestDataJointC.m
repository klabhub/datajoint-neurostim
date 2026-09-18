classdef TestDataJointC < TestDataJointPipelineBase
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
    end
end
