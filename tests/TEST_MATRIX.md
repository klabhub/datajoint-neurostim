# Regression test matrix

This file records the intended coverage of the repository's regression tests.
Update it when tests are added, removed, or their scope changes. It describes
behavioral intent; generated code-coverage reports provide execution details.

| Component or path | Test file | Coverage type | Main behaviors covered | Known gaps |
| --- | --- | --- | --- | --- |
| `nsScan` | `TestScan.m` | Dependency-free unit test | Empty/nonexistent dates, Neurostim filename parsing, date-folder scanning, filters, byte metadata | JSON metadata, database insertion, file-content parsing |
| `nsScan` JSON metadata | `TestScan.m` + `fixtures/nsScanJson` | Dependency-free fixture test | Subject, session, and experiment definitions and metadata, lock flags, analyze filtering | JSON validation failures, alternate metadata tags |
| `readJson` | `TestUtilities.m` | Unit test | Scalar number conversion and empty files | Nested/invalid JSON |
| `catstruct` | `TestUtilities.m` | Unit test | Compatible defaults for missing fields | Nested struct edge cases |
| `resampleTrials` | `TestUtilities.m` | Unit test | Condition-preserving selection and partitioning | Boundary and weighted cases |
| `retimeWithNan` | `TestUtilities.m` | Unit test | Strict and partial-missing sum behavior | Other retiming methods |
| `nsInitializeDataJoint` | `TestDataJointPipeline.m` | Integration test | Temporary project/schema initialization | Invalid names and multiple packages |
| `ns.File` | `TestDataJointPipeline.m` | DataJoint integration test | File population from a synthetic experiment | Missing/excluded files and checksum errors |
| `ns.C` / `ns.CChannel` | `TestDataJointPipeline.m` | DataJoint integration test | Synthetic reader contract, signal/time storage, channels, sampling rate | Reader failures and malformed outputs |
| `ns.Dimension` / parts | `TestDataJointPipeline.m` | DataJoint integration test | Plugin-driven conditions and trial expansion | Multiple dimensions and restrictions |
| `ns.Epoch` / `ns.EpochChannel` | `TestDataJointPipeline.m` | DataJoint integration test | Alignment, epoch window, channel rows, extracted signal, key source | Artifact/plugin attrition edge cases |
| `chunkedDelete` | `TestDataJointPipeline.m` | DataJoint integration test | Batched deletion of EpochChannel rows and preservation of the master Epoch row | External FK discovery, reconnect/retry path, confirmation prompt |
| `ns.Tepoch` / `ns.TepochChannel` | `TestDataJointPipeline.m` | DataJoint integration test | FFT transformation with `average={}`, per-trial/per-channel output, dependent/independent metadata, and unaveraged row counts | Channel/trial/condition averaging, restrictions, alternate transforms, and window edge cases |

## Running tests

Run the normal suite with:

```matlab
addpath("tests", "-begin")
addpath(pwd)
addpath(fullfile(pwd, "datajoint-matlab"))
addpath(fullfile(pwd, "tools", "mym", mexext))
results = run_tests;
assertSuccess(results)
```

The DataJoint integration tests create and remove a uniquely named database.
They require a DataJoint account with permission to create and drop databases.

## Coverage reports

Use `run_coverage` to generate a MATLAB HTML report in `coverage/html`. The report is generated
output and should not be committed; update this matrix when test intent changes.
