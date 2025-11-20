# Package index

## Main Functionality

Fit the model with
[`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md),
visualize the results with
[`plot.SLOPE()`](https://jolars.github.io/SLOPE/reference/plot.SLOPE.md),
and produce predictions for new data with
[`predict.SLOPE()`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md).

- [`SLOPE()`](https://jolars.github.io/SLOPE/reference/SLOPE.md) :
  Sorted L-One Penalized Estimation
- [`plot(`*`<SLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/plot.SLOPE.md)
  : Plot Coefficients
- [`print(`*`<SLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/print.SLOPE.md)
  [`print(`*`<TrainedSLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/print.SLOPE.md)
  : Print Results from SLOPE Fit
- [`summary(`*`<SLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/summary.SLOPE.md)
  : Summarize SLOPE Model
- [`print(`*`<summary_SLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/print.summary_SLOPE.md)
  : Print Summary of SLOPE Model
- [`predict(`*`<SLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md)
  [`predict(`*`<GaussianSLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md)
  [`predict(`*`<BinomialSLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md)
  [`predict(`*`<PoissonSLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md)
  [`predict(`*`<MultinomialSLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/predict.SLOPE.md)
  : Generate Predictions from SLOPE Models
- [`coef(`*`<SLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/coef.SLOPE.md)
  : Obtain Coefficients
- [`deviance(`*`<SLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/deviance.SLOPE.md)
  : Model Deviance
- [`score()`](https://jolars.github.io/SLOPE/reference/score.md) :
  Compute One of Several Loss Metrics on a New Data Set

## Clusters

Analyze the cluster strucure from SLOPE fits

- [`plotClusters()`](https://jolars.github.io/SLOPE/reference/plotClusters.md)
  : Plot Cluster Structure

## Model Tuning

Use [`cvSLOPE()`](https://jolars.github.io/SLOPE/reference/cvSLOPE.md)
to perform cross-validation for selecting the optimal regularization
parameters,
[`plot.TrainedSLOPE()`](https://jolars.github.io/SLOPE/reference/plot.TrainedSLOPE.md)
to visualize the results, and
[`refit()`](https://jolars.github.io/SLOPE/reference/refit.md) to fit
the final model to your data set.

- [`cvSLOPE()`](https://jolars.github.io/SLOPE/reference/cvSLOPE.md) :
  Tune SLOPE with Cross-Validation
- [`trainSLOPE()`](https://jolars.github.io/SLOPE/reference/trainSLOPE.md)
  : Train a SLOPE Model
- [`plot(`*`<TrainedSLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/plot.TrainedSLOPE.md)
  : Plot Results from Cross-Validation
- [`print(`*`<SLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/print.SLOPE.md)
  [`print(`*`<TrainedSLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/print.SLOPE.md)
  : Print Results from SLOPE Fit
- [`refit()`](https://jolars.github.io/SLOPE/reference/refit.md) : Refit
  SLOPE Model with Optimal Parameters
- [`summary(`*`<TrainedSLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/summary.TrainedSLOPE.md)
  : Summarize TrainedSLOPE Model
- [`print(`*`<summary_TrainedSLOPE>`*`)`](https://jolars.github.io/SLOPE/reference/print.summary_TrainedSLOPE.md)
  : Print Summary of TrainedSLOPE Model

## Utilities

Helper functions for various tasks.

- [`plotDiagnostics()`](https://jolars.github.io/SLOPE/reference/plotDiagnostics.md)
  : Plot Results from Diagnostics Collected During Model Fitting
- [`regularizationWeights()`](https://jolars.github.io/SLOPE/reference/regularizationWeights.md)
  : Generate Regularization (Penalty) Weights for SLOPE
- [`sortedL1Prox()`](https://jolars.github.io/SLOPE/reference/sortedL1Prox.md)
  : Sorted L1 Proximal Operator

## Data Sets

Data sets buddled with the package.

- [`abalone`](https://jolars.github.io/SLOPE/reference/abalone.md) :
  Abalone
- [`bodyfat`](https://jolars.github.io/SLOPE/reference/bodyfat.md) :
  Bodyfat
- [`heart`](https://jolars.github.io/SLOPE/reference/heart.md) : Heart
  Disease
- [`student`](https://jolars.github.io/SLOPE/reference/student.md) :
  Student Performance
- [`wine`](https://jolars.github.io/SLOPE/reference/wine.md) : Wine
  Cultivars
- [`glioma`](https://jolars.github.io/SLOPE/reference/glioma.md) :
  Glioma Metabolomics
