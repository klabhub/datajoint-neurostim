classdef TestDataJointDeletePlan < TestDataJointPipelineBase
    % These tests use only the suite's isolated temporary database.
    methods (Test)
        function stopsAtEmptyChildOfNonemptyRoot(testCase)
            insert(ns.DimensionParm,struct('dimension','deletePlanEmpty', ...
                'paradigm','synthetic','parms',struct( ...
                'plg','synthetic','prm','condition','atTrialTime',0)));
            root = ns.DimensionParm & 'dimension="deletePlanEmpty"';

            [rels,counts,capped,approximate] = dj.planDelete(root);

            testCase.verifyEqual(counts(1),1);
            testCase.verifyEqual(counts(2:end),zeros(size(counts(2:end))));
            testCase.verifyFalse(any(capped | approximate));
            testCase.verifyTrue(all(arrayfun(@(rel) contains(rel.sql,'FALSE'),rels(3:end))));
        end

        function retainsNonemptyBranches(testCase)
            root = ns.Dimension & testCase.experimentKey();

            [rels,counts,capped,approximate] = dj.planDelete(root);

            names = string({rels.className});
            testCase.verifyEqual(counts(names=="ns.Dimension"),1);
            testCase.verifyEqual(counts(names=="ns.DimensionCondition"),2);
            testCase.verifyEqual(counts(names=="ns.DimensionTrial"),3);
            testCase.verifyFalse(any(capped | approximate));
        end

        function emptyRootPrunesEntireCascade(testCase)
            root = ns.DimensionParm & 'dimension="deletePlanMissing"';

            [rels,counts,capped,approximate] = dj.planDelete(root);

            testCase.verifyEqual(counts,zeros(size(counts)));
            testCase.verifyFalse(any(capped | approximate));
            testCase.verifyTrue(all(arrayfun(@(rel) contains(rel.sql,'FALSE'),rels(2:end))));
        end

        function largeRootRetainsCappedCount(testCase)
            root = ns.PluginParameter & testCase.experimentKey();

            [~,counts,capped,approximate] = dj.planDelete(root,1);

            testCase.verifyEqual(counts(1),1);
            testCase.verifyTrue(capped(1));
            testCase.verifyFalse(approximate(1));
        end
    end
end
