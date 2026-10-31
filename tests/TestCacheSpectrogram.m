classdef TestCacheSpectrogram < matlab.unittest.TestCase
    properties (TestParameter)
        configuration = struct( ...
            'shortEpoch',struct('n',188,'fs',250,'limits',[0.5 50],'voices',6,'boundary',"reflection",'precision',"double"), ...
            'longEpoch',struct('n',1875,'fs',250,'limits',[0.5 50],'voices',6,'boundary',"reflection",'precision',"double"), ...
            'samplingRate',struct('n',512,'fs',100,'limits',[2 35],'voices',6,'boundary',"reflection",'precision',"double"), ...
            'frequencyRange',struct('n',512,'fs',250,'limits',[5 20],'voices',6,'boundary',"reflection",'precision',"double"), ...
            'voices',struct('n',512,'fs',250,'limits',[5 20],'voices',12,'boundary',"reflection",'precision',"double"), ...
            'periodic',struct('n',512,'fs',250,'limits',[5 20],'voices',6,'boundary',"periodic",'precision',"double"), ...
            'singlePrecision',struct('n',512,'fs',250,'limits',[5 20],'voices',6,'boundary',"reflection",'precision',"single"))
        transform = struct( ...
            'cwt',struct('name',"spectrogram",'options',struct('limits',[0.05 0.4],'voicesPerOctave',4)), ...
            'cohen',struct('name',"wavelet",'options',struct('limits',[0.05 0.4],'nfrex',4,'fwhm',[20 3])))
    end

    methods (TestClassSetup)
        function requireWaveletToolbox(testCase)
            testCase.assumeTrue(~isempty(which('cwtfilterbank')) && ...
                license('test','Wavelet_Toolbox'));
            testsFolder = fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(testsFolder)));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(testsFolder));
        end
    end

    methods (Test)
        function agreesWithNativeCwt(testCase,configuration)
            c = configuration;
            time = (0:c.n-1)'/c.fs;
            signal = cast(sin(2*pi*10*time)+0.3*cos(2*pi*23*time),c.precision);
            fb = cwtfilterbank(SignalLength=c.n,SamplingFrequency=c.fs, ...
                Wavelet='amor',FrequencyLimits=c.limits, ...
                VoicesPerOctave=c.voices,Boundary=c.boundary);
            [expected,frequency] = wt(fb,signal);

            [actual,idv] = ns.cache.do_spectrogram(signal,c.fs, ...
                limits=c.limits,voicesPerOctave=c.voices,boundary=c.boundary);

            testCase.verifyEqual(idv,"xt");
            testCase.verifyEqual(actual.Properties.VariableNames,{'power','xt'});
            testCase.verifyEqual(actual.power{1},flipud(abs(expected).^2).', ...
                'AbsTol',1e-6);
            testCase.verifyEqual(actual.xt{1},flipud(frequency(:)).', ...
                'AbsTol',1e-12);
            testCase.verifyEqual(actual.xt{2},time.','AbsTol',1e-12);
            testCase.verifySize(actual.power{1},[c.n numel(frequency)]);
            testCase.verifyClass(actual.power{1},char(c.precision));
        end

        function repeatedCallsTransformNewSignal(testCase)
            signal = sin(2*pi*10*(0:1023)'/250);
            first = ns.cache.do_spectrogram(signal,250,limits=[5 20]);

            second = ns.cache.do_spectrogram(2*signal,250,limits=[5 20]);

            testCase.verifyEqual(second.power{1},4*first.power{1},'AbsTol',1e-10);
            testCase.verifyEqual(second.xt,first.xt);
        end

        function powerUsesSquaredCwtMagnitude(testCase)
            signal = 2*cos(2*pi*10*(0:2499)'/250);

            actual = ns.cache.do_spectrogram(signal,250,limits=[5 10]);

            testCase.verifyEqual(actual.xt{1}(end),10,'AbsTol',1e-12);
            testCase.verifyEqual(mean(actual.power{1}(751:1750,end)),4, ...
                'AbsTol',0.01);
        end

        function reportsActualShortEpochFrequencyRange(testCase)
            signal = sin(2*pi*10*(0:187)'/250);

            actual = ns.cache.do_spectrogram(signal,250,limits=[0.5 50]);

            testCase.verifyGreaterThan(min(actual.xt{1}),0.5);
            testCase.verifyLessThanOrEqual(max(actual.xt{1}),50);
            testCase.verifyGreaterThan(diff(actual.xt{1}),0);
            testCase.verifySize(actual.power{1},[188 numel(actual.xt{1})]);
        end

        function rejectsTooShortSignal(testCase)
            testCase.verifyError(@() ns.cache.do_spectrogram(ones(3,1),250), ...
                'ns:cache:SpectrogramSignalLength');
        end

        function rejectsReversedFrequencyLimits(testCase)
            testCase.verifyError(@() ns.cache.do_spectrogram(ones(128,1),250, ...
                limits=[20 5]),'ns:cache:SpectrogramFrequencyLimits');
        end

        function rejectsFrequenciesAboveNyquist(testCase)
            testCase.verifyError(@() ns.cache.do_spectrogram(ones(128,1),250, ...
                limits=[5 126]),'ns:cache:SpectrogramFrequencyLimits');
        end

        function normalizesEachFrequencyToBaseline(testCase)
            time = (0:1023)/250-2;
            signal = (1+0.2*time').*sin(2*pi*10*time');
            raw = ns.cache.do_spectrogram(signal,250,time,limits=[5 20]);
            baseline = time >= -1.5 & time <= -0.5;
            expected = raw.power{1}./mean(raw.power{1}(baseline,:),1);

            actual = ns.cache.do_spectrogram(signal,250,time, ...
                limits=[5 20],baselineWindow=[-1500 -500]);
            scaled = ns.cache.do_spectrogram(2*signal,250,time, ...
                limits=[5 20],baselineWindow=[-1500 -500]);

            testCase.verifyEqual(actual.power{1},expected,'AbsTol',1e-12);
            testCase.verifyEqual(mean(actual.power{1}(baseline,:),1), ...
                ones(1,size(expected,2)),'AbsTol',1e-12);
            testCase.verifyEqual(scaled.power,actual.power,'AbsTol',1e-12);
            testCase.verifyEqual(actual.xt,raw.xt);
        end

        function recomputesBaselineForEachWindow(testCase)
            time = (0:1023)/250;
            signal = (1+time').*sin(2*pi*10*time');
            ns.cache.do_spectrogram(signal,250,limits=[5 20],baselineWindow=[0 500]);
            raw = ns.cache.do_spectrogram(signal,250,limits=[5 20]);

            actual = ns.cache.do_spectrogram(signal,250, ...
                limits=[5 20],baselineWindow=[2000 3000]);

            baseline = time >= 2 & time <= 3;
            testCase.verifyEqual(actual.power{1}, ...
                raw.power{1}./mean(raw.power{1}(baseline,:),1),'AbsTol',1e-12);
        end

        function zeroBaselineProducesNaNs(testCase)
            actual = ns.cache.do_spectrogram(zeros(512,1),250, ...
                limits=[5 20],baselineWindow=[0 500]);

            testCase.verifyTrue(all(isnan(actual.power{1}),'all'));
        end

        function rejectsUnavailableBaseline(testCase)
            testCase.verifyError(@() ns.cache.do_spectrogram(ones(512,1),250, ...
                baselineWindow=[-500 -1]),'ns:cache:SpectrogramBaselineWindow');
            testCase.verifyError(@() ns.cache.do_spectrogram(ones(512,1),250, ...
                baselineWindow=[500 0]),'ns:cache:SpectrogramBaselineWindow');
            testCase.verifyError(@() ns.cache.do_spectrogram(ones(512,1),250, ...
                baselineWindow=[0 100 500]),'ns:cache:SpectrogramBaselineWindow');
        end

        function computeUsesEpochBaselineAfterCropping(testCase,transform)
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('NS_PARPOOL','0'));
            time = -128:127;
            signal = (2+time/256).*sin(2*pi*0.1*time);
            cache = testsupport.VariableLengthCache(table({signal},{[-128 127 256]}, ...
                'VariableNames',{'signal','time'}));
            options = transform.options;
            options.baselineWindow = [-90000 -10000];

            actual = cache.compute(struct(transform.name,options), ...
                average=string.empty,timeWindow=[-100000 100000]);

            selected = time >= -100 & time <= 100;
            transformFunction = str2func("ns.cache.do_" + transform.name);
            rawOptions = namedargs2cell(transform.options);
            raw = transformFunction(signal(selected)',1,time(selected),rawOptions{:});
            baseline = raw.xt{2} >= -90 & raw.xt{2} <= -10;
            expected = raw.power{1}./mean(raw.power{1}(baseline,:),1);
            testCase.verifyEqual(actual.(transform.name+"_power"){1},expected,'AbsTol',1e-12);
            testCase.verifyEqual(actual.(transform.name+"_xt"),raw.xt);
        end

        function baselineNormalizationMatchesForBothTransforms(testCase,transform)
            time = -128:127;
            signal = (2+time'/256).*sin(2*pi*0.1*time');
            transformFunction = str2func("ns.cache.do_"+transform.name);
            options = namedargs2cell(transform.options);
            raw = transformFunction(signal,1,time,options{:});
            baseline = time >= -90 & time <= -10;

            actual = transformFunction(signal,1,time,options{:},baselineWindow=[-90000 -10000]);
            scaled = transformFunction(2*signal,1,time,options{:},baselineWindow=[-90000 -10000]);

            testCase.verifyEqual(actual.power{1}, ...
                raw.power{1}./mean(raw.power{1}(baseline,:),1),'AbsTol',1e-12);
            testCase.verifyEqual(mean(actual.power{1}(baseline,:),1), ...
                ones(1,size(raw.power{1},2)),'AbsTol',1e-12);
            testCase.verifyEqual(scaled.power,actual.power,'AbsTol',1e-12);
            testCase.verifyEqual(actual.xt,raw.xt);
        end

        function waveletHandlesZeroAndMissingBaseline(testCase)
            actual = ns.cache.do_wavelet(zeros(256,1),1, ...
                limits=[0.05 0.4],fwhm=[20 3],nfrex=4,baselineWindow=[0 50000]);

            testCase.verifyTrue(all(isnan(actual.power{1}),'all'));
            testCase.verifyError(@() ns.cache.do_wavelet(ones(256,1),1, ...
                baselineWindow=[-500 -1]),'ns:cache:WaveletBaselineWindow');
            testCase.verifyError(@() ns.cache.do_wavelet(ones(256,1),1, ...
                baselineWindow=[500 0]),'ns:cache:WaveletBaselineWindow');
        end

        function computePreservesTimeFrequencyContract(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('NS_PARPOOL','0'));
            signal = sin(2*pi*0.1*(0:255));
            cache = testsupport.VariableLengthCache(table({signal},{[0 255 256]}, ...
                'VariableNames',{'signal','time'}));
            options = struct('limits',[0.05 0.4],'voicesPerOctave',4);

            [actual,mapping] = cache.compute(struct('spectrogram',options), ...
                average=string.empty);
            expected = ns.cache.do_spectrogram(signal',1, ...
                limits=options.limits,voicesPerOctave=options.voicesPerOctave);

            testCase.verifyEqual(actual.spectrogram_power,expected.power);
            testCase.verifyEqual(actual.spectrogram_xt,expected.xt);
            testCase.verifyEqual(mapping("spectrogram_xt"),"spectrogram_power");
            testCase.verifyEqual(cache.independent,"time");
            testCase.verifyEqual(cache.dependent,"signal");
        end

        function computeHandlesMultipleGroups(testCase,transform)
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('NS_PARPOOL','0'));
            signal = sin(2*pi*0.1*(0:255));
            condition = repelem((1:12)',4);
            trial = repmat([1;1;2;2],12,1);
            channel = repmat([1;2;1;2],12,1);
            data = table(ones(48,1),condition,trial,channel, ...
                repmat({signal},48,1),repmat({[0 255 256]},48,1), ...
                'VariableNames',{'subject','condition','trial','channel','signal','time'});
            cache = testsupport.VariableLengthCache(data);
            referenceCache = testsupport.VariableLengthCache(table({signal},{[0 255 256]}, ...
                'VariableNames',{'signal','time'}));
            fun = struct(transform.name,transform.options);
            reference = referenceCache.compute(fun,average=string.empty);

            [actual,mapping] = cache.compute(fun,average=["trial" "channel"]);

            powerName = transform.name + "_power";
            axesName = transform.name + "_xt";
            testCase.verifyEqual(height(actual),12);
            testCase.verifyEqual(actual.(powerName),repmat(reference.(powerName),12,1));
            testCase.verifyEqual(actual.(axesName),repmat(reference.(axesName),12,1));
            testCase.verifyEqual(mapping(axesName),powerName);
            testCase.verifyEqual(actual.nrtrials,2*ones(12,1));
            testCase.verifyEqual(actual.nrchannels,2*ones(12,1));
        end
    end
end
