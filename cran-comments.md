## Test environments

* Ubuntu 20.04, R-release (local)
* Windows, R-devel (win-builder)
* Ubuntu Linux 16.04 LTS, R-release, GCC (rhub)
* Fedora Linux, R-devel, clang, gfortran (rhub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (rhub) 
* OSX R-devel (github)

## R CMD check results

0 errors | 0 warnings | 1 notes

> Found the following (possibly) invalid URLs:
>   URL: http://www.jstor.org/stable/2346178
>     From: inst/doc/introduction.html
>     Status: 403
>     Message: Forbidden

I believe the status message above is a false positive; the URL works as it should.

> Days since last update: 5

This a second patch (see below).

## Resubmission

This is a resubmission. In this version I have actually fixed the UBSAN errors (<https://cran.r-project.org/web/checks/check_results_SLOPE.html>). I can no longer reproduce them using the ASAN/UBSAN docker from rhub/rocker (<https://builder.r-hub.io/status/original/SLOPE_0.3.1.9000.tar.gz-15e2ec2c51bd4992bcd4060f5727bb6f>).

I have also fixed the test that failed for the windows build test check.
