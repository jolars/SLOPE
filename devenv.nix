{
  pkgs,
  lib,
  config,
  inputs,
  ...
}:

{
  packages = with pkgs; [
    bashInteractive
    autoconf
    go-task
    quartoMinimal
    pandoc
    llvmPackages.openmp
  ];

  languages = {
    r = {
      enable = true;

      package = (
        pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            BH
            Matrix
            Rcpp
            RcppEigen
            SparseM
            bigmemory
            caret
            covr
            devtools
            dplyr
            e1071
            glmnet
            here
            knitr
            languageserver
            rhub
            rmarkdown
            scales
            spelling
            testthat
            tidyverse
            usethis
          ];
        }
      );
    };
  };
}
