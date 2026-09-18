classdef TestPreprocessTiming < matlab.unittest.TestCase
    properties (TestParameter)
        sampleRate = {250,512,2000,30000}
    end

    methods (Test)
        function resamplingPreservesDurationAndTone(testCase,sampleRate)
            time = (0:2*sampleRate-1)'/sampleRate;
            signal = sin(2*pi*10*time);
            targetRate = sampleRate/2;

            [actual,actualTime] = prep.preprocess(signal,time, ...
                struct('resample',struct('frequency',targetRate)));

            testCase.verifySize(actual,[sampleRate 1]);
            testCase.verifyEqual(actualTime,(0:sampleRate-1)'/targetRate,AbsTol=1e-10);
            % Avoid filter transients at the edges of the recording.
            interior = actualTime > 0.2 & actualTime < 1.8;
            testCase.verifyEqual(actual(interior),sin(2*pi*10*actualTime(interior)), ...
                AbsTol=0.01);
        end

        function decimationPreservesSampleSpacingAndOrigin(testCase)
            time = 2 + (0:999)'/1000;
            signal = sin(2*pi*10*time);

            [actual,actualTime] = prep.preprocess(signal,time, ...
                struct('decimate',struct('frequency',250)));

            testCase.verifySize(actual,[250 1]);
            testCase.verifyEqual(actualTime,2+(0:249)'/250,AbsTol=1e-10);
        end
    end
end
