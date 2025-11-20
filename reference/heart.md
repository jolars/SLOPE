# Heart Disease

Diagnostic attributes of patients classified as having heart disease or
not.

## Usage

``` r
heart
```

## Format

270 observations from 17 variables represented as a list consisting of a
binary factor response vector `y`, with levels 'absence' and 'presence'
indicating the absence or presence of heart disease and `x`: a sparse
feature matrix of class 'dgCMatrix' with the following variables:

- age:

  age

- bp:

  diastolic blood pressure

- chol:

  serum cholesterol in mg/dl

- hr:

  maximum heart rate achieved

- old_peak:

  ST depression induced by exercise relative to rest

- vessels:

  the number of major blood vessels (0 to 3) that were colored by
  fluoroscopy

- sex:

  sex of the participant: 0 for male, 1 for female

- angina:

  a dummy variable indicating whether the person suffered
  angina-pectoris during exercise

- glucose_high:

  indicates a fasting blood sugar over 120 mg/dl

- cp_typical:

  typical angina

- cp_atypical:

  atypical angina

- cp_nonanginal:

  non-anginal pain

- ecg_abnormal:

  indicates a ST-T wave abnormality (T wave inversions and/or ST
  elevation or depression of \> 0.05 mV)

- ecg_estes:

  probable or definite left ventricular hypertrophy by Estes' criteria

- slope_flat:

  a flat ST curve during peak exercise

- slope_downsloping:

  a downwards-sloping ST curve during peak exercise

- thal_reversible:

  reversible defect

- thal_fixed:

  fixed defect

## Source

Dua, D. and Karra Taniskidou, E. (2017). UCI Machine Learning Repository
<http://archive.ics.uci.edu/ml/>. Irvine, CA: University of California,
School of Information and Computer Science.

<https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html#heart>

## Preprocessing

The original dataset contained 13 variables. The nominal of these were
dummycoded, removing the first category. No precise information
regarding variables `chest_pain`, `thal` and `ecg` could be found, which
explains their obscure definitions here.

## See also

Other datasets:
[`abalone`](https://jolars.github.io/SLOPE/reference/abalone.md),
[`bodyfat`](https://jolars.github.io/SLOPE/reference/bodyfat.md),
[`student`](https://jolars.github.io/SLOPE/reference/student.md),
[`wine`](https://jolars.github.io/SLOPE/reference/wine.md)
