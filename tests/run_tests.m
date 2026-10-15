function results = run_tests
%RUN_TESTS Run the regression tests for the main repository.

testsFolder = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(testsFolder);
addpath(repositoryRoot);

suite = matlab.unittest.TestSuite.fromFolder(testsFolder, ...
    'IncludingSubfolders', true);
results = run(suite);

if nargout == 0
    assertSuccess(results);
end
end
