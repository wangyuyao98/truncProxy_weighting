## Test environments

- Local macOS (Apple Silicon), R 4.5.1
- GitHub Actions:
  - ubuntu-latest, R-release
  - windows-latest, R-release
  - macos-latest, R-release

## R CMD check results

- `R CMD check --as-cran truncProxy_0.1.0.tar.gz`

Local `R CMD check --as-cran` completed with 3 NOTEs and no errors or warnings.

The NOTEs were:

- This is a new submission.
- `README.md` and `NEWS.md` could not be checked locally because `pandoc` is not installed.
- HTML validation could not be fully checked locally because `HTML Tidy` is not recent enough.

These are local environment issues and not package issues.

GitHub Actions checks passed on Linux, macOS, and Windows.

## Downstream dependencies

- This is a new package, so there are no downstream dependencies.
