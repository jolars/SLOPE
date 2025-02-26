{
  description = "A basic flake with a shell";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  inputs.systems.url = "github:nix-systems/default";
  inputs.flake-utils = {
    url = "github:numtide/flake-utils";
    inputs.systems.follows = "systems";
  };

  outputs =
    { nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      {
        devShells.default = pkgs.mkShell {
          shellHook = ''
            mkdir -p "$(pwd)/_libs"
            export R_LIBS_USER="$(pwd)/_libs"
          '';
          packages =
            let
              SLOPE = (
                pkgs.rPackages.buildRPackage {
                  name = "SLOPE";
                  src = ./.;
                  propagatedBuildInputs = with pkgs.rPackages; [
                    foreach
                    ggplot2
                    Matrix
                    Rcpp
                    RcppEigen
                    bench
                    covr
                    dplyr
                    knitr
                    rmarkdown
                    scales
                    spelling
                    testthat
                    tidyr
                    SparseM
                    caret
                    e1071
                  ];
                }
              );
            in
            with pkgs;
            [
              bashInteractive
              autoconf
              go-task
              (rWrapper.override {
                packages = with rPackages; [
                  devtools
                  languageserver
                  SLOPE
                  tidyverse
                  usethis
                  rhub
                ];
              })
            ];
        };
      }
    );
}
