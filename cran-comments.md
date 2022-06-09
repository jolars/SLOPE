## Test environments

* Fedora 36, R-release, gcc (local)
* Windows, R-devel (win-builder)
* Windows, R-release (win-builder)
* Debian Linux, R-devel, GCC ASAN/UBSAN (rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (rhub)
* Ubuntu 20.04, R-release, R-devel (github)
* Windows Server 2022, R-release (github)
* Mac OS X 11.6.6, R-release (github)

## R CMD check results

0 errors | 0 warnings | 1 notes

>> * checking CRAN incoming feasibility ... NOTE
>> Maintainer: 'Johan Larsson <johan.larsson@stat.lu.se>'
>> Found the following (possibly) invalid URLs:
>>   URL: https://doi.org/10.1111/j.1541-0420.2007.00843.x
>>     From: man/SLOPE.Rd
>>     Status: 503
>>     Message: Service Unavailable
>>   URL: https://doi.org/10.1214/15-AOAS842
>>     From: inst/doc/introduction.html
>>     Status: 500
>>     Message: Internal Server Error
>>   URL: https://doi.org/10/gfgwzt
>>     From: man/SLOPE.Rd
>>     Status: 500
>>     Message: Internal Server Error
>>   URL: https://www.jstor.org/stable/2346178
>>     From: inst/doc/introduction.html
>>     Status: 403
>>     Message: Forbidden
>> 
>> Found the following (possibly) invalid DOIs:
>>   DOI: 10/gfgwzt
>>     From: DESCRIPTION
>>           inst/CITATION
>>     Status: Internal Server Error
>>     Message: 500

I believe these notes to be false positives. All of the links work when I
check them.
