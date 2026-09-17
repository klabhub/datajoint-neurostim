classdef TestDataJointSnr < TestDataJointPipelineBase
    methods (Test)
        function fullSnrExampleFindsInjectedPeakNearTenHz(testCase)
            originalPool = getenv('NS_PARPOOL');
            testCase.addTeardown(@() setenv('NS_PARPOOL',originalPool));
            setenv('NS_PARPOOL','0');
            originalRng = rng;
            testCase.addTeardown(@() rng(originalRng));
            rng(42,'twister');

            % A long epoch resolves the offset from the 10 Hz search target.
            fs = 250;
            injectedFrequency = 9.8;
            t = (0:30*fs-1)'/fs;
            signal = sin(2*pi*injectedFrequency*t) + 0.05*randn(size(t));
            etag = 'snrExampleEpoch';
            insert(ns.EpochParm,struct('etag',etag,'ctag','synthetic', ...
                'dimension','condition','window',1000*[t(1) t(end)], ...
                'align',struct('plugin','synthetic','event','startTime')));
            testCase.addTeardown(@() TestDataJointSnr.removeEpoch(etag));
            key = fetch(proj(ns.C)*proj(ns.Dimension) & testCase.experimentKey());
            testCase.assertNumElements(key,1);
            key.etag = etag;
            insert(ns.Epoch,mergestruct(key,struct( ...
                'time',[1000*t(1) 1000*t(end) numel(t)], ...
                'prep',struct(),'art',struct(),'plg',struct())));
            insert(ns.EpochChannel,mergestruct(key,struct( ...
                'channel',1,'trial',1,'onset',0,'signal',signal)));
            epochs = ns.EpochChannel & key;
            testCase.assertEqual(count(epochs),1);

            fun.pspectrum = struct('FrequencyLimits',[0 50]);
            fun.snr = struct('signalHalfWidth',1,'noiseHalfWidth',2);
            fun.peak = struct('searchFrequencies',10 , ...
                'searchRangeHalfWidth',1);
            [result,dependent,independent] = compute(epochs,fun,average=string.empty);

            testCase.verifyEqual(dependent,["peakFrequency" "magnitude"]);
            testCase.verifyEqual(independent,"searchFrequency");
            testCase.assertEqual(height(result),1);
            searchFrequency = cell2mat(result.searchFrequency);
            peakFrequency = cell2mat(result.peakFrequency);
            magnitude = cell2mat(result.magnitude);
            testCase.verifyEqual(searchFrequency(:),[2;6;10;24]);
            testCase.assertNumElements(peakFrequency,4);
            testCase.assertNumElements(magnitude,4);
            atTenHz = searchFrequency == 10;
            testCase.verifyEqual(peakFrequency(atTenHz),injectedFrequency,AbsTol=0.1);
            testCase.verifyTrue(isfinite(magnitude(atTenHz)));
            testCase.verifyGreaterThan(magnitude(atTenHz),100);
        end
    end

    methods (Static, Access=private)
        function removeEpoch(etag)
            previousSafeMode = dj.config('safemode');
            restoreSafeMode = onCleanup(@() dj.config('safemode',previousSafeMode));
            dj.config('safemode',false);
            delete(ns.EpochParm & struct('etag',etag));
        end
    end
end
