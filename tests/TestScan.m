classdef TestScan < matlab.unittest.TestCase
    % Dependency-light tests for nsScan file discovery behavior.

    properties
        temporaryRoot
    end

    methods (TestMethodSetup)
        function createTemporaryRoot(testCase)
            testCase.temporaryRoot = string(tempname);
            mkdir(testCase.temporaryRoot);
        end
    end

    methods (TestMethodTeardown)
        function removeTemporaryRoot(testCase)
            if ~isempty(testCase.temporaryRoot) && isfolder(testCase.temporaryRoot)
                rmdir(testCase.temporaryRoot, 's');
            end
        end
    end

    methods (Test)
        function emptyDayReturnsEmptyTables(testCase)
            [subjects, sessions, experiments, locked] = nsScan( ...
                'root', testCase.temporaryRoot, ...
                'date', datetime(2024, 1, 15), ...
                'schedule', "d", ...
                'verbose', false);

            testCase.verifyEmpty(subjects);
            testCase.verifyEmpty(sessions);
            testCase.verifyEmpty(experiments);
            testCase.verifyEqual(locked, struct( ...
                'experiment', false, ...
                'session', false, ...
                'subject', false));
        end

        function nonexistentDayDoesNotCreateResults(testCase)
            [~, ~, experiments] = nsScan( ...
                'root', testCase.temporaryRoot, ...
                'date', datetime(2030, 12, 31), ...
                'schedule', "d", ...
                'verbose', false);

            testCase.verifyEmpty(experiments);
        end

        function scansNeurostimFilesFromDateFoldersWithoutReadingContents(testCase)
            sessionDate = datetime(2024, 2, 29);
            dayFolder = fullfile(testCase.temporaryRoot, '2024', '02', '29');
            testCase.assertTrue(mkdir(dayFolder));

            fileNames = [
                "1.behavior.090000.mat"
                "1.stimulation.101530.mat"
                "2.behavior.111500.mat"
            ];
            for fileName = fileNames'
                fileID = fopen(fullfile(dayFolder, fileName), 'w');
                testCase.assertNotEqual(fileID, -1);
                fclose(fileID);
            end

            [subjects, sessions, experiments] = nsScan( ...
                'root', testCase.temporaryRoot, ...
                'date', sessionDate, ...
                'schedule', "d", ...
                'paradigm', ["behavior" "stimulation"], ...
                'readJson', false, ...
                'readFileContents', false, ...
                'verbose', false);

            testCase.verifyEqual(subjects.subject, ["1"; "2"]);
            testCase.verifyEqual(sessions.session_date, repmat("2024-02-29", 2, 1));
            testCase.verifyEqual(sessions.subject, ["1"; "2"]);
            testCase.verifyEqual(experiments.session_date, repmat("2024-02-29", 3, 1));
            testCase.verifyEqual(experiments.file, fileNames);
            testCase.verifyEqual(experiments.starttime, ["09:00:00"; "10:15:30"; "11:15:00"]);
            testCase.verifyEqual(experiments.subject, ["1"; "1"; "2"]);
            testCase.verifyEqual(experiments.paradigm, ["behavior"; "stimulation"; "behavior"]);
            testCase.verifyEqual(experiments.bytes, zeros(3, 1));
        end
    end
end
