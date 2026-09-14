classdef TestUtilities < matlab.unittest.TestCase
    % Regression tests for utilities in the main repository.

    properties
        temporaryFile
    end

    methods (TestMethodTeardown)
        function removeTemporaryFile(testCase)
            if ~isempty(testCase.temporaryFile) && isfile(testCase.temporaryFile)
                delete(testCase.temporaryFile);
            end
        end
    end

    methods (Test)
        function readJsonConvertsScalarNumbersToStrings(testCase)
            testCase.temporaryFile = [tempname '.json'];
            writeJson(struct('number', 42, 'text', 'hello'), ...
                testCase.temporaryFile);

            actual = readJson(testCase.temporaryFile);

            testCase.verifyEqual(actual.number, "42");
            testCase.verifyEqual(actual.text, "hello");
        end

        function readJsonReturnsEmptyStructForEmptyFile(testCase)
            testCase.temporaryFile = [tempname '.json'];
            fclose(fopen(testCase.temporaryFile, 'w'));

            actual = readJson(testCase.temporaryFile);

            testCase.verifyEmpty(actual);
        end

        function catstructAddsMissingFieldsWithCompatibleDefaults(testCase)
            first = struct('number', 1, 'name', "first");
            second = struct('number', 2, 'extra', 3.5);

            actual = catstruct(2, first, second);

            testCase.verifyEqual([actual.number], [1 2]);
            testCase.verifyEqual(actual(1).name, "first");
            testCase.verifyEqual(actual(2).name, "");
            testCase.verifyEqual(actual(1).extra, []);
            testCase.verifyEqual(actual(2).extra, 3.5);
        end

        function resampleTrialsPreservesConditionsAndPartitionsTrials(testCase)
            rng(7);
            condition = [1 1 1 2 2 2 3 3 3];

            [selected, omitted] = resampleTrials(condition, false, 0.5);

            testCase.verifyEqual(unique(condition(selected)), [1 2 3]);
            testCase.verifyEmpty(intersect(selected, omitted));
            testCase.verifyEqual(sort([selected omitted]), 1:numel(condition));
        end

        function retimeWithNanStrictlyPropagatesMissingValuesForSum(testCase)
            times = seconds([0 0 1])';
            input = timetable(times, [1; NaN; 3], ...
                'VariableNames', {'value'});

            output = retimeWithNan(input, seconds(0:2), 'sum', ...
                'AnyNanIsNanSum', true);

            testCase.verifyTrue(isnan(output.value(1)));
        end

        function retimeWithNanCanIgnorePartialMissingValuesForSum(testCase)
            times = seconds([0 0 1])';
            input = timetable(times, [1; NaN; 3], ...
                'VariableNames', {'value'});

            output = retimeWithNan(input, seconds(0:2), 'sum', ...
                'AnyNanIsNanSum', false);

            testCase.verifyEqual(output.value(1), 1);
        end

        function fftAcceptsNamedOptions(testCase)
            signal = {sin(2*pi*(0:4)'/5)};
            result = ns.cache.do_fft(signal,1000,n=8);

            testCase.verifyEqual(result.frequency{1},(0:4)*1000/8);
            testCase.verifySize(result.amplitude{1},[5 1]);
        end

        function pspectrumAcceptsNamedOptions(testCase)
            signal = {sin(2*pi*(0:7)'/8)};
            result = ns.cache.do_pspectrum(signal,1000,FrequencyLimits=[0 250], ...
                Leakage=1,TwoSided=false);

            testCase.verifyEqual(result.frequency{1}(1),0);
            testCase.verifyLessThanOrEqual(result.frequency{1}(end),250);
            testCase.verifySize(result.power{1},size(result.frequency{1}));
        end
        function pmtmAcceptsStructOptions(testCase)
            signal = {sin(2*pi*(0:31)'/8)};
            result = ns.cache.do_pmtm(signal,nw=4,nfft=16,fs=100);

            testCase.verifySize(result.frequency{1},[9 1]);
            testCase.verifySize(result.power{1},size(result.frequency{1}));
        end

        function waveletAcceptsStructOptions(testCase)
            signal = {sin(2*pi*(0:31)'/8)};
            result = ns.cache.do_wavelet(signal,100, ...
                nfrex=4,fwhm=[2 1],limits=[5 20]);

            testCase.verifyEqual(result.frequency{1},[5 10 15 20]');
            testCase.verifySize(result.power{1},[4 32]);
            testCase.verifySize(result.time{1},[1 32]);
        end

        function snrAcceptsStructOptions(testCase)
            signal = {1 + (0:16)'/16};
            frequencies = (0:16)'/4;
            result = ns.cache.do_snr(signal,frequencies, ...
                signalHalfWidth=1,noiseHalfWidth=2);

            testCase.verifyEqual(result.frequency{1},frequencies);
            testCase.verifySize(result.snr{1},[17 1]);
        end

        function searchPeaksAcceptsStructOptions(testCase)
            signal = [(0:7)' (7:-1:0)'];
            frequencies = (0:7)';
            result = ns.cache.do_search_peaks(signal,frequencies, ...
                searchFrequencies=[2 6],searchRangeHalfWidth=0.5);

            testCase.verifyEqual(result.searchFrequency{1},[2 2; 6 6]);
            testCase.verifyEqual(result.peakFrequency{1},[2 2; 6 6]);
            testCase.verifyEqual(result.magnitude{1},[2 5; 6 1]);
        end    end
end
