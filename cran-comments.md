## Test environments

- Fedora 40 gcc, R-release (local)
- Windows, R-devel (win-builder)
- Windows, R-release (win-builder)
- Ubuntu 22.04, R-release, R-devel (github)
- Windows Server 2022 10.0.20348, R-release (github)
- Mac OS X 14.5, R-release (github)
- Ubuntu 22.04.4 with rchk, R-devel (rhub)
- Fedora Linux 38 with valgrind, R-devel (rhub)
- Ubuntu 22.04.4 LTS with clang-asan, R-devel (rhub)

## R CMD check results

0 errors | 0 warnings | 1 notes

* checking CRAN incoming feasibility ... [17s] NOTE
  Maintainer: 'Johan Larsson <johanlarsson@outlook.com>'

  New maintainer:
    Johan Larsson <johanlarsson@outlook.com>
  Old maintainer(s):
    Johan Larsson <johan.larsson@stat.lu.se>

## Reverse Dependencies

We checked reverse dependencies geneSLOPE and sgs for compatibility with the new
version and found no issues.
