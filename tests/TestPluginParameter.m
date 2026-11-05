classdef TestPluginParameter < matlab.unittest.TestCase
    properties (TestParameter)
        blobPayload = struct('scalar',42,'rowVector',[1 2 3], ...
            'columnVector',[1;2;3],'matrix',reshape(1:6,2,3), ...
            'array3d',reshape(1:8,2,2,2),'text','example', ...
            'cellValue',{{1,'text'}},'empty',[])
    end
    methods (Test)
        function boundaryPreservesEachBlobPayload(testCase,blobPayload)
            tuple = struct('plugin_name','sample','property_name','a', ...
                'property_type','Event','property_value',{blobPayload}, ...
                'property_time',{blobPayload},'property_trial',{blobPayload}, ...
                'property_nstime',{blobPayload});
            tuples = repmat(tuple,2,1);
            tuples(2).property_name = 'b';

            actual = ns.PluginParameter.parameterTable(tuples);

            testCase.verifyClass(actual.property_value,'cell');
            testCase.verifyClass(actual.property_time,'cell');
            testCase.verifyClass(actual.property_trial,'cell');
            testCase.verifyClass(actual.property_nstime,'cell');
            testCase.verifyEqual(actual.property_value,{blobPayload;blobPayload});
            testCase.verifyEqual(actual.property_time,{blobPayload;blobPayload});
            testCase.verifyEqual(actual.property_trial,{blobPayload;blobPayload});
            testCase.verifyEqual(actual.property_nstime,{blobPayload;blobPayload});
            testCase.verifyEqual(actual.property_name,["a";"b"]);
            testCase.verifyEqual(actual.Properties.VariableNames,fieldnames(tuples)');
        end

        function unpacksScalarBlobRows(testCase)
            names = ["firstFrame";"trialStopTime"];

            actual = ns.PluginParameter.nestedPluginTable(names,{100;500}, ...
                {-Inf;-Inf},{100;500},{1;1},["Global";"Global"]);

            testCase.verifyEqual(actual.firstframe.data,100);
            testCase.verifyEqual(actual.trialstoptime.data,500);
            testCase.verifyEqual(actual.trialstoptime.clocktime,500);
            testCase.verifyEqual(actual.trialstoptime.trial,1);
        end

        function preservesMatrixRowPayloads(testCase)
            actual = ns.PluginParameter.nestedPluginTable(["a";"b"], ...
                {[1 2 3];[4 5 6]},{[0 1 2];[3 4 5]}, ...
                {[10 11 12];[13 14 15]},{[1 1 1];[2 2 2]},["Event";"Event"]);

            testCase.verifyEqual(actual.a.data,[1 2 3]);
            testCase.verifyEqual(actual.b.data,[4 5 6]);
            testCase.verifyEqual(actual.b.trialtime,[3 4 5]);
            testCase.verifyEqual(actual.b.clocktime,[13 14 15]);
        end

        function preservesSinglePropertyVector(testCase)
            actual = ns.PluginParameter.nestedPluginTable("a", ...
                {[1;2;3]},{[0;1;2]},{[10;11;12]},{[1;1;1]},"Event");

            testCase.verifyEqual(actual.a.data,[1;2;3]);
            testCase.verifyEqual(actual.a.clocktime,[10;11;12]);
        end

        function preservesCellPayloads(testCase)
            actual = ns.PluginParameter.nestedPluginTable(["a";"b"], ...
                {[1 2];[3 4 5]},{[0 1];[2 3 4]}, ...
                {[10 11];[12 13 14]},{[1 1];[2 2 2]},["Event";"Event"]);

            testCase.verifyEqual(actual.a.data,[1 2]);
            testCase.verifyEqual(actual.b.data,[3 4 5]);
        end

        function cellValuedPropertyProducesScalarStruct(testCase)
            payload = {42,[1 2 3]};

            actual = ns.PluginParameter.nestedPluginTable("a", ...
                {payload},{[0 1]},{[100 101]},{[1 1]},"Event");

            testCase.verifySize(actual.a,[1 1]);
            testCase.verifyEqual(actual.a.data,payload);
        end

        function decodesByteStreamWithoutExtraCellLayer(testCase)
            payload = {42,[1 2 3]};
            encoded = getByteStreamFromArray(payload);

            actual = ns.PluginParameter.nestedPluginTable("a", ...
                {encoded},{[0 1]},{[100 101]},{[1 1]},"ByteStream");

            testCase.verifySize(actual.a,[1 1]);
            testCase.verifyEqual(actual.a.data,payload);
        end

        function rejectsInputsOutsideBlobContract(testCase)
            testCase.verifyError(@() ns.PluginParameter.nestedPluginTable( ...
                "a",42,{0},{100},{1},"Global"),'ns:PluginParameter:BlobShape');
        end
    end
end
