# Student performance

A data set of the attributes of 382 students in secondary education
collected from two schools. The goal is to predict the grade in math and
Portugese at the end of the third period. See the cited sources for
additional information.

## Usage

``` r
student
```

## Format

382 observations from 13 variables represented as a list consisting of a
binary factor response matrix `y` with two responses: `portugese` and
`math` for the final scores in period three for the respective subjects.
The list also contains `x`: a sparse feature matrix of class 'dgCMatrix'
with the following variables:

- school_ms:

  student's primary school, 1 for Mousinho da Silveira and 0 for Gabriel
  Pereira

- sex:

  sex of student, 1 for male

- age:

  age of student

- urban:

  urban (1) or rural (0) home address

- large_family:

  whether the family size is larger than 3

- cohabitation:

  whether parents live together

- Medu:

  mother's level of education (ordered)

- Fedu:

  fathers's level of education (ordered)

- Mjob_health:

  whether the mother was employed in health care

- Mjob_other:

  whether the mother was employed as something other than the specified
  job roles

- Mjob_services:

  whether the mother was employed in the service sector

- Mjob_teacher:

  whether the mother was employed as a teacher

- Fjob_health:

  whether the father was employed in health care

- Fjob_other:

  whether the father was employed as something other than the specified
  job roles

- Fjob_services:

  whether the father was employed in the service sector

- Fjob_teacher:

  whether the father was employed as a teacher

- reason_home:

  school chosen for being close to home

- reason_other:

  school chosen for another reason

- reason_rep:

  school chosen for its reputation

- nursery:

  whether the student attended nursery school

- internet:

  Pwhether the student has internet access at home

## Source

P. Cortez and A. Silva. Using Data Mining to Predict Secondary School
Student Performance. In A. Brito and J. Teixeira Eds., Proceedings of
5th FUture BUsiness TEChnology Conference (FUBUTEC 2008) pp. 5-12,
Porto, Portugal, April, 2008, EUROSIS, ISBN 978-9077381-39-7.

Dua, D. and Karra Taniskidou, E. (2017). UCI Machine Learning Repository
<http://archive.ics.uci.edu/ml/>. Irvine, CA: University of California,
School of Information and Computer Science.

## Preprocessing

All of the grade-specific predictors were dropped from the data set.
(Note that it is not clear from the source why some of these predictors
are specific to each grade, such as which parent is the student's
guardian.) The categorical variables were dummy-coded. Only the final
grades (G3) were kept as dependent variables, whilst the first and
second period grades were dropped.

## See also

Other datasets:
[`abalone`](https://jolars.github.io/SLOPE/reference/abalone.md),
[`bodyfat`](https://jolars.github.io/SLOPE/reference/bodyfat.md),
[`heart`](https://jolars.github.io/SLOPE/reference/heart.md),
[`wine`](https://jolars.github.io/SLOPE/reference/wine.md)
