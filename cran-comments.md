## Submission

This is a new submission.

## Test environments

* Linux (x86_64, conda-forge gcc 16.2), R 4.4.2: R CMD check --as-cran
  (with --run-donttest, PDF/HTML manual, remote URL checks)
* GitHub Actions (R-CMD-check.yaml): macOS, Windows, Ubuntu R-devel/release/oldrel-1
  -- to be run after pushing
* win-builder (R-devel, R-release) and macOS builder -- to be run before submission

## R CMD check results

0 errors | 0 warnings | 3 notes (local Linux)

* checking CRAN incoming feasibility ... NOTE
  New submission. The README R-CMD-check badge URL returned 404 because the
  GitHub Actions workflow had not yet been pushed; it resolves once the
  workflow has run.
* checking for future file timestamps ... NOTE ("unable to verify current
  time") -- a property of the local sandbox, not the package.
* checking compilation flags used ... NOTE ('-march=nocona') -- set by the
  local conda-forge toolchain, not by the package (src/Makevars only sets
  PKG_LIBS = -pthread).

## Notes for CRAN

* Examples and tests use at most 2 threads.
* Functions that contact NCBI web services are wrapped in \donttest{} and
  their tests are skipped on CRAN.
* BLAST+ is optional (SystemRequirements); the BLAST engine is only used if
  the user asks for it, and its tests are skipped when BLAST+ is absent.
