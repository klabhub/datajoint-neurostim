classdef TestCacheInsert < matlab.unittest.TestCase
    properties (TestParameter)
        validSignal = {1:5,42,ones(3,4),ones(2,3,4)}
        invalidSignal = {(1:5)',[],zeros(0,3),{1:5}}
    end
    methods (Test)
        function acceptsValidPayload(testCase,validSignal)
            target = testsupport.CacheInsertValidator;
            testCase.verifyWarningFree(@() insert(target,struct('signal',validSignal)));
        end
        function rejectsInvalidPayload(testCase,invalidSignal)
            target = testsupport.CacheInsertValidator;
            tuple.signal = invalidSignal;
            testCase.verifyError(@() insert(target,tuple),'ns:cache:InvalidSignalShape');
        end
        function checksEveryTuple(testCase)
            target = testsupport.CacheInsertValidator;
            tuples = struct('signal',{1:3,(1:3)'});
            testCase.verifyError(@() insert(target,tuples),'ns:cache:InvalidSignalShape');
        end
        function validatesCellInput(testCase)
            target = testsupport.CacheInsertValidator;
            testCase.verifyWarningFree(@() insert(target,{1:3;4:6}));
            testCase.verifyError(@() insert(target,{1:3;(4:6)'}), ...
                'ns:cache:InvalidSignalShape');
        end
        function acceptsEmptyBatch(testCase)
            target = testsupport.CacheInsertValidator;
            testCase.verifyWarningFree(@() insert(target,struct('signal',{})));
        end
    end
end
