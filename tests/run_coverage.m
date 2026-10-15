function results = run_coverage
%RUN_COVERAGE Run regression tests and write a MATLAB HTML coverage report.

% Reset the DataJoint MEX connection and persistent package/schema functions so repeated runs use fresh test databases.
try
    dj.conn.close;
catch
    % No connection has been established yet.
end
clear dj.conn
clear classes
clear functions

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
