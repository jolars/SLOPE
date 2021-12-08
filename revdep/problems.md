# grpSLOPE

<details>

* Version: 0.2.1
* Source code: https://github.com/cran/grpSLOPE
* URL: https://github.com/agisga/grpSLOPE.git
* BugReports: https://github.com/agisga/grpSLOPE/issues
* Date/Publication: 2016-11-20 09:18:04
* Number of recursive dependencies: 40

Run `revdep_details(,"grpSLOPE")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    Running examples in ‘grpSLOPE-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: coef.grpSLOPE
    > ### Title: Extract model coefficients
    > ### Aliases: coef.grpSLOPE
    > 
    > ### ** Examples
    > 
    > set.seed(1)
    > A   <- matrix(rnorm(100^2), 100, 100)
    > grp <- rep(rep(letters[1:20]), each=5)
    > b   <- c(rep(1, 20), rep(0, 80))
    > y   <- A %*% b + rnorm(10) 
    > result <- grpSLOPE(X=A, y=y, group=grp, fdr=0.1)
    Warning: 'SLOPE::prox_sorted_L1' is deprecated.
    See help("Deprecated")
    Error in sorted_l1_prox(x, lambda) : Not a matrix.
    Calls: grpSLOPE ... proximalGradientSolverGroupSLOPE -> proxGroupSortedL1 -> <Anonymous> -> sorted_l1_prox
    Execution halted
    ```

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      ══ testthat results  ══════════════════════════════════════════════════════════════════════════════════════════════════════════════
      [ OK: 66 | SKIPPED: 0 | WARNINGS: 25 | FAILED: 25 ]
      1. Error: (unknown) (@test_generics.R#11) 
      2. Error: when the groups are consecutive blocks (@test_grpSLOPE.R#22) 
      3. Error: when the groups are not consecutive blocks (@test_grpSLOPE.R#39) 
      4. Error: with non-zero intercept (@test_grpSLOPE.R#56) 
      5. Error: when the groups are consecutive blocks (@test_grpSLOPE.R#140) 
      6. Error: when the groups are not consecutive blocks (@test_grpSLOPE.R#158) 
      7. Error: with non-zero intercept (@test_grpSLOPE.R#177) 
      8. Error: when rank of group submatrix is smaller than the group size (@test_grpSLOPE.R#198) 
      9. Error: when the groups are consecutive blocks (@test_grpSLOPE.R#232) 
      1. ...
      
      Error: testthat unit tests failed
      Execution halted
    ```

# owl

<details>

* Version: 0.1.1
* Source code: https://github.com/cran/owl
* URL: https://github.com/jolars/owl, https://jolars.github.io/owl
* BugReports: https://github.com/jolars/owl/issues
* Date/Publication: 2020-02-11 10:50:08 UTC
* Number of recursive dependencies: 105

Run `revdep_details(,"owl")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      5/5 mismatches (average diff: 2.09)
      [1] 2.33 - 0.233 == 2.09
      [2] 2.33 - 0.233 == 2.09
      [3] 2.33 - 0.233 == 2.09
      [4] 2.33 - 0.233 == 2.09
      [5] 2.33 - 0.233 == 2.09
      
      ══ testthat results  ══════════════════════════════════════════════════════════════════════════════════════════════════════════════
      [ OK: 55 | SKIPPED: 0 | WARNINGS: 4 | FAILED: 3 ]
      1. Failure: SLOPE and owl agree for gaussian designs (@test-slope.R#21) 
      2. Failure: SLOPE and owl agree when computing lambda sequences (@test-slope.R#37) 
      3. Failure: SLOPE and owl agree when computing lambda sequences (@test-slope.R#37) 
      
      Error: testthat unit tests failed
      Execution halted
    ```

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 11.4Mb
      sub-directories of 1Mb or more:
        libs  10.9Mb
    ```

