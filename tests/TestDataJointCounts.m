classdef TestDataJointCounts < TestDataJointPipelineBase
    methods (Test)
        function shortPreviewProvidesExactCount(testCase)
            rel = testCase.displayRelation(1000,12);
            [text,queries] = TestDataJointCounts.captureDisplay(rel);
            testCase.verifyTrue(contains(text,'2 tuples'));
            testCase.verifyFalse(contains(text,'>2 tuples'));
            testCase.verifyEqual(queries,1);
        end

        function emptyPreviewProvidesExactZero(testCase)
            rel = testCase.displayRelation(1000,12) & 'channel=-1';
            [text,queries] = TestDataJointCounts.captureDisplay(rel);
            testCase.verifyTrue(contains(text,'0 tuples'));
            testCase.verifyEqual(queries,1);
        end

        function previewCanProveCountCapExceeded(testCase)
            rel = testCase.displayRelation(1,12);
            [text,queries] = TestDataJointCounts.captureDisplay(rel);
            testCase.verifyTrue(contains(text,'>1 tuples'));
            testCase.verifyEqual(queries,1);
        end

        function fullPreviewStillCountsWhenNecessary(testCase)
            rel = testCase.displayRelation(1000,1);
            [text,queries] = TestDataJointCounts.captureDisplay(rel);
            testCase.verifyTrue(contains(text,'2 tuples'));
            testCase.verifyEqual(queries,2);
        end

        function boundedCountDistinguishesExactThreshold(testCase)
            rel = testCase.displayRelation(1000,12);
            [below,cappedBelow] = countUpTo(rel,1);
            [equal,cappedEqual] = countUpTo(rel,2);
            [above,cappedAbove] = countUpTo(rel,3);
            [unlimited,cappedUnlimited] = countUpTo(rel,Inf);
            testCase.verifyEqual([below equal above unlimited],[1 2 2 2]);
            testCase.verifyEqual([cappedBelow cappedEqual cappedAbove cappedUnlimited], ...
                [true false false false]);
        end
    end
    methods (Access=private)
        function rel = displayRelation(testCase,countLimit,displayLimit)
            previousCount = dj.config('countLimit');
            previousDisplay = dj.config('displayLimit');
            previousEnabled = dj.config('displayCount');
            testCase.addTeardown(@() dj.config('countLimit',previousCount));
            testCase.addTeardown(@() dj.config('displayLimit',previousDisplay));
            testCase.addTeardown(@() dj.config('displayCount',previousEnabled));
            dj.config('countLimit',countLimit);
            dj.config('displayLimit',displayLimit);
            dj.config('displayCount',true);
            rel = ns.CChannel & testCase.experimentKey();
            testCase.assertNotEmpty(rel.header.primaryKey); % Warm metadata caches.
        end
    end
    methods (Static, Access=private)
        function [text,queries] = captureDisplay(rel)
            before = rel.schema.conn.query('SHOW SESSION STATUS LIKE ''Questions''');
            text = evalc('disp(rel)');
            after = rel.schema.conn.query('SHOW SESSION STATUS LIKE ''Questions''');
            queries = str2double(after.Value{1})-str2double(before.Value{1})-1;
        end
    end
end
