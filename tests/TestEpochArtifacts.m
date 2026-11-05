classdef TestEpochArtifacts < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addRepository(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function masksOverlappingWindowsWithPadding(testCase)
            signal = repmat(sin(2*pi*5*(0:299)'/100),[1 2 2]);
            signal(125,1,:) = 1000;
            options = TestEpochArtifacts.options();
            options.epoch_no = [7;20];

            [actual,badBy,details] = prep.epochArtifacts(signal,100,options);

            testCase.verifySize(actual,size(signal));
            testCase.verifyTrue(all(isnan(actual(41:210,1,:)),'all'));
            testCase.verifyEqual(actual(1:40,1,:),signal(1:40,1,:));
            testCase.verifyEqual(actual(211:end,1,:),signal(211:end,1,:));
            testCase.verifyEqual(actual(:,2,:),signal(:,2,:));
            testCase.verifyEmpty(badBy.all);
            testCase.verifyEqual(details.intervals{:,:},[7 41 210]);
            testCase.verifyEqual(details.maskedFraction,[170/300;0],'AbsTol',1e-12);
            testCase.verifyEqual(details.windowFlags.badByAmplitude,int32([2;3]));
        end

        function includesFinalPartialStride(testCase)
            signal = sin(2*pi*5*(0:324)'/100);
            signal(end) = 1000;
            options = TestEpochArtifacts.options();
            options.padding = 0;

            [actual,~,details] = prep.epochArtifacts(signal,100,options);

            testCase.verifyTrue(all(isnan(actual(226:325))));
            testCase.verifyEqual(actual(1:225),signal(1:225));
            testCase.verifyEqual(details.intervals{:,:},[1 226 325]);
        end

        function missingWindowsDoNotEnterReferenceStatistics(testCase)
            signal = sin(2*pi*5*(0:299)'/100);
            signal(125) = NaN;
            options = TestEpochArtifacts.options();
            options.padding = 0;

            [actual,~,details] = prep.epochArtifacts(signal,100,options);

            testCase.verifyEqual(details.windowFlags.badByMissing,int32([2;3]));
            testCase.verifyTrue(all(isnan(actual(51:200))));
            testCase.verifyEqual(actual(201:end),signal(201:end));
        end

        function handlesShortAndEntirelyMissingEpochs(testCase)
            options = TestEpochArtifacts.options();

            [actual,badBy,details] = prep.epochArtifacts(NaN(30,1),100,options);

            testCase.verifyTrue(all(isnan(actual)));
            testCase.verifyEmpty(badBy.all);
            testCase.verifyEqual(details.intervals{:,:},[1 1 30]);
            testCase.verifyEqual(details.maskedFraction,1);
        end

        function honorsOriginalTrialExclusions(testCase)
            signal = repmat(sin(2*pi*5*(0:299)'/100),1,2);
            signal(125,:) = 1000;
            options = TestEpochArtifacts.options();
            options.epoch_no = [7;20];
            options.exclude = 7;

            [actual,~,details] = prep.epochArtifacts(signal,100,options);

            testCase.verifyEqual(actual(:,1),signal(:,1));
            testCase.verifyEqual(details.maskedFraction(1),0);
            testCase.verifyEqual(unique(details.windows.trial),20);
            testCase.verifyTrue(any(isnan(actual(:,2))));
        end

        function disabledDetectionPreservesSignal(testCase)
            signal = ones(100,1);
            signal(20) = NaN;
            options = TestEpochArtifacts.options();
            options.enable = false;

            [actual,badBy,details] = prep.epochArtifacts(signal,100,options);

            testCase.verifyEqual(actual,signal);
            testCase.verifyEmpty(badBy.all);
            testCase.verifyEmpty(details.intervals);
            testCase.verifyEqual(details.maskedFraction,0);
        end

        function defaultActionMatchesExistingDetector(testCase)
            signal = repmat(sin(2*pi*5*(0:299)'/100),1,2);
            signal(125,2) = 1000;
            options = rmfield(TestEpochArtifacts.options(), ...
                {'action','window','step','padding'});
            options.epoch_no = [7;20];
            args = namedargs2cell(options);
            expected = prep.artifactDetection(permute(signal,[2 3 1]),100,args{:});

            [actual,badBy,details] = prep.epochArtifacts(signal,100,options);

            testCase.verifyEqual(actual,signal);
            testCase.verifyEqual(badBy.all,expected.all);
            testCase.verifyEqual(badBy.badByAmplitude,int32(20));
            testCase.verifyEqual(details.action,"reject");
        end

        function rejectsUncoveredSamples(testCase)
            options = TestEpochArtifacts.options();
            options.step = 2000;

            testCase.verifyError(@() prep.epochArtifacts(ones(300,1),100,options), ...
                'prep:epochArtifacts:Step');
        end
    end

    methods (Static, Access=private)
        function options = options()
            options = struct('enable',true,'action',"nan",'window',1000, ...
                'step',500,'padding',100,'flat_threshold_sd',0, ...
                'amplitude_threshold_peak',50,'amplitude_channelfrac',0.5, ...
                'variance_z_threshold',Inf,'hf_cutoff_hz',20, ...
                'hf_z_threshold',Inf,'correlation_z_threshold',Inf);
        end
    end
end
