classdef TestEeglabPreprocess < matlab.unittest.TestCase
    % Regression tests for ephys.eeglab.preprocess.

    methods (Test)
        function numberedStepsAreAppliedInOrder(testCase)
            repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
            mockFolder = fullfile(repositoryRoot, 'tests', 'fixtures', 'eeglab');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mockFolder, ...
                'IncludingSubfolders', false));
            clear pop_reref eeg_checkset

            EEG = struct('etc', struct());
            parms.eeglab.reref1 = true;
            parms.eeglab.reref2 = true;

            actual = ephys.eeglab.preprocess(EEG, parms);

            testCase.verifyEqual(actual.etc.rerefCalls, [1 2]);
        end

        function fullRankIcaReplacesInfinitePcaWithRank(testCase)
            repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
            mockFolder = fullfile(repositoryRoot, 'tests', 'fixtures', 'eeglab');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(mockFolder, ...
                'IncludingSubfolders', false));
            clear pop_runica eeg_checkset

            EEG = struct('data', eye(3, 4), 'pnts', 4, 'etc', struct());
            parms.eeglab.ica = struct('pca', inf);

            actual = ephys.eeglab.preprocess(EEG, parms);

            testCase.verifyEqual(actual.etc.runicaArgs, {'pca', 3});
        end
    end
end
