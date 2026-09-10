function results = run_coverage
%RUN_COVERAGE Run regression tests and write a MATLAB HTML coverage report.

testsFolder = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(testsFolder);
addpath(repositoryRoot);
addpath(fullfile(repositoryRoot,'datajoint-matlab'));
addpath(fullfile(repositoryRoot,'tools','mym',mexext));

coverageFolder = fullfile(testsFolder,'coverage','html');
if ~isfolder(coverageFolder)
    mkdir(coverageFolder);
end

suite = matlab.unittest.TestSuite.fromFolder(testsFolder, ...
    'IncludingSubfolders',true);
runner = matlab.unittest.TestRunner.withTextOutput;
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport
runner.addPlugin(CodeCoveragePlugin.forFolder(repositoryRoot, ...
    'Producing',CoverageReport(coverageFolder)));

results = runner.run(suite);
if nargout == 0
    assertSuccess(results);
end
end
