## Test environments

* local ubuntu 19.10, R 3.6.2
* win-builder (devel and release)
* Debian Linux, R-devel, GCC ASAN/UBSAN (on rhub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (on rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (on rhub)
* Ubuntu Linux 16.04 LTS, R-devel with rchk (on rhub)
* windows (release), OSX (release, devel), and linux (release) on github

## R CMD check results

0 errors | 0 warnings | 1 note

> * checking CRAN incoming feasibility ... NOTE
> Maintainer: 'Johan Larsson <johan.larsson@stat.lu.se>'
> 
> New maintainer:
>   Johan Larsson <johan.larsson@stat.lu.se>
> Old maintainer(s):
>   Evan Patterson <epatters@stanford.edu>

Evan Patterson has transferred maintainership to me and has
sent the CRAN team an e-mail notifying you of this change, which I quote here:

> Transfer of maintainership of SLOPE package
> 
> Dear CRAN admins,
> 
> I am writing to authorize transfer of maintainership of the SLOPE package 
> from me to Malgorzata Bogdan and/or Johan Larsson.
> 
> https://cran.r-project.org/web/packages/SLOPE/index.html
> 
> Thank you,
> Evan

> * checking CRAN incoming feasibility ... NOTE
> Possibly mis-spelled words in DESCRIPTION:
>   Bogdan (31:6)
>   al (31:16)
>   et (31:13)
  
This is a false positive. These words are correctly spelled.

> * checking CRAN incoming feasibility ... NOTE
> Found the following (possibly) invalid URLs:
>   URL: http://www.jstor.org/stable/2346178
>     From: inst/doc/introduction.html
>     Status: 403
>     Message: Forbidden
>   URL: https://doi.org/10.1137/080716542
>     From: man/SLOPE.Rd
>     Status: Error
>     Message: libcurl error code 56:
>       	Recv failure: Connection was reset

I believe these are also false positives. I can access both of these
URLs without problem.

  
### Resubmission

This is a resubmission. In this version I have modifed some examples
in order to make them run faster and pass the automated checks.
