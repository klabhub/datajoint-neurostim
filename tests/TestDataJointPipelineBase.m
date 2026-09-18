classdef (Abstract, SharedTestFixtures = {TestDataJointPipelineFixture}) TestDataJointPipelineBase < matlab.unittest.TestCase
    % Shared integration-test fixture for the ns pipeline.
    properties
        projectRoot
        dataRoot
        databaseName
    end
    properties (Constant, Access = private)
        subject = 'joe'
        sessionDate = '2024-02-29'
        startTime = '09:00:00'
    end
    methods (TestClassSetup)
        function useSharedDatabase(testCase)
            fixture = testCase.getSharedTestFixtures();
            testCase.projectRoot = fixture.projectRoot;
            testCase.dataRoot = fixture.dataRoot;
            testCase.databaseName = fixture.databaseName;
        end
    end
    methods (Access = protected)
        function report(~,message,varargin)
            fprintf('[TestDataJointPipeline] %s\n',sprintf(message,varargin{:}));
        end
        function key = experimentKey(testCase)
            key = struct('subject',testCase.subject, ...
                'session_date',testCase.sessionDate, ...
                'starttime',testCase.startTime);
        end
    end
end
