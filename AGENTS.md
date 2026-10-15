# AGENTS.md

## Project conventions

- This is a MATLAB project using DataJoint, with package code primarily under `+ns`, `+ephys`, `+prep`, `+sbx`, and `+fc`.
- Preserve existing DataJoint table declarations, foreign-key relationships, defaults, and schema comments unless the task explicitly requires schema changes.
- Follow the existing MATLAB class and package structure. Prefer small, focused changes over broad refactors.
- Preserve existing uncommitted user changes. Never reset, discard, or overwrite unrelated modifications.
- Use `apply_patch` for source edits.

## Testing

- The regression suite is in `tests/run_tests.m`.
- Run it with:

  ```matlab
  addpath("tests", "-begin")
  addpath(pwd)
  addpath(fullfile(pwd, "datajoint-matlab"))
  addpath(fullfile(pwd, "tools", "mym", mexext))
  results = run_tests;
  assertSuccess(results)

##  DataJoint safety
- Do not run destructive database commands manually unless explicitly requested.
- Prefer the test suite’s isolated temporary database for integration testing.
- Do not modify production schemas or data as part of ordinary tests.
- Verify cleanup after interrupted integration tests when practical.

## Communication
- For implementation requests, make the change and run the most relevant tests.
- For diagnostic requests, explain the cause before changing code.
- Keep updates concise and evidence-based.
- Distinguish test failures from environment or infrastructure failures.
