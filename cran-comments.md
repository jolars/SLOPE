## Test environments

* Ubuntu 20.04, R-release (local)
* Windows, R-devel (win-builder)
* Ubuntu Linux 16.04 LTS, R-release, GCC (rhub)
* Fedora Linux, R-devel, clang, gfortran (rhub)
* OSX R-devel (github)

## R CMD check results

0 errors | 0 warnings | 1 notes

> Found the following (possibly) invalid URLs:
>   URL: http://www.jstor.org/stable/2346178
>     From: inst/doc/introduction.html
>     Status: 403
>     Message: Forbidden

I believe the status message above is a false positive; the URL works as it should.

> Days since last update: 4

This is a patch to fix the broken install on solaris as well as rectify ASAN/UBSAN errors (<https://cran.r-project.org/web/checks/check_results_SLOPE.html>).
