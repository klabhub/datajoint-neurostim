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
    end
end
