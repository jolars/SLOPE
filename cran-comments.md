## Test environments

- nixos-unstable gcc, R-release (local)
- Windows, R-devel (win-builder)
- Windows, R-release (win-builder)
- Ubuntu-latest, R-release and devel, R-devel (github)
- Windows-latest, R-release (github)
- Mac OS X-latest, R-release (github)
- Fedora 38, valgrind, R-devel (rhub)
- Uuntu 22.04.5, clang-asan, R-devel (rhub)
- Uuntu 22.04.5, clang-ubsan, R-devel (rhub)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse Dependencies

We checked reverse dependencies geneSLOPE and sgs for compatibility with the new
version and found no issues.

## Bug fix for M1mac system

This release fixes a test failure on the M1mac system
<https://www.stats.ox.ac.uk/pub/bdr/M1mac/SLOPE.out>. Testing on mac-builder
<<https://mac.R-project.org/macbuilder/results/1751374240-4c650b290d8cc31c/> and
rhub <https://github.com/jolars/SLOPE/actions/runs/15998323131> suggests that
the issue is now resolved.
