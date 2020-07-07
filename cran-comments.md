## Test environments

* ubuntu 20.04, r-release (local)
* windows R-devel (on win-builder)
* ubuntu linux 16.04 LTS, r-release, GCC (on rhub)
* Fedora Linux, R-devel, clang, gfortran (on rhub)
* OSX r-devel (on github)


## R CMD check results

0 errors | 0 warnings | 1 notes

winbuilder is giving me the following note but
I think that it is a false positive. The URL works fine.

> Found the following (possibly) invalid URLs:
>   URL: http://www.jstor.org/stable/2346178
>     From: inst/doc/introduction.html
>     Status: 403
>     Message: Forbidden

## Reverse dependencies

SLOPE has no reverse dependencies.
